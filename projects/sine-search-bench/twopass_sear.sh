#!/usr/bin/env bash
# Two-pass SINE search with sear's geometry and filters:
#   pass 1  fasta36 (ktup 6) on 2.5 kb / 500 bp sliding windows, parts of <=20k windows  -> candidate loci
#   pass 2a ssearch36 (sear's flags, -Z fixed) on candidate regions (+-PAD bp)            -> full SW alignments
#   pass 2b ssearch36 (sear's flags) on the genome with candidates masked (optional)      -> rescue of pass-1 misses
# Output: <out>/loci.bed (chrom start end strand pctid bits source), after sear's 65% / 0.8-length filter.
# Usage: twopass_sear.sh genome.fa query.fa outdir [THREADS=4] [PAD=400] [RESIDUAL=1]
set -euo pipefail
G=$1; Q=$2; O=$3; T=${4:-4}; PAD=${5:-400}; RES=${6:-1}; MODE=${MODE:-regions}; RULE=${RULE:-hsp}; GAP=${GAP:-50}; RES_WIN=${RES_WIN:-2500}; RES_OVL=${RES_OVL:-500}; RES_Z=${RES_Z:-11}
SS=${SSEARCH36:-/home/user/src/fasta36/bin/ssearch36}; FA=${FASTA36:-/home/user/src/fasta36/bin/fasta36}
B=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$O"; L=$(grep -v '^>' "$Q" | tr -d '\n' | wc -c)
t0=$(date +%s)
python3 - "$G" "$O/win.fa" <<'PY'
import sys; sys.path.insert(0,sys.argv[0] if False else '/home/user/knowlege/projects/sine-search-bench')
from bench import chunk_fasta, split_parts
chunk_fasta(sys.argv[1], sys.argv[2], 2500, 500); print(len(split_parts(sys.argv[2], 20000)), "parts")
PY
t1=$(date +%s)
for p in "$O"/win.fa.part*.fa; do $FA -g -3 -T $T -Q -n -z 11 -E 2 -m 8C "$Q" "$p" 6; done > "$O/pass1.raw" 2>/dev/null
t2=$(date +%s); echo "pass1 fasta36: $((t2-t1)) s, $(grep -vc '^#' "$O/pass1.raw") HSPs"
python3 - "$O/pass1.raw" "$G" "$O" "$PAD" <<'PY'
import sys, collections; sys.path.insert(0,'/home/user/knowlege/projects/sine-search-bench')
from simulate import read_fasta
raw, G, O, pad = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
hits = collections.defaultdict(list)
for l in open(raw):
    if l.startswith('#') or not l.strip(): continue
    f = l.rstrip('\n').split('\t'); c, rng = f[1].rsplit(':', 1); off = int(rng.split('-')[0])
    s, e = sorted((int(f[8]), int(f[9]))); hits[c].append((s - 1 + off, e + off))
n = tot = 0
with open(f"{O}/regions.fa", "w") as o, open(f"{O}/regions.bed", "w") as b:
    for c, seq in read_fasta(G):
        hs = sorted(hits.get(c, [])); merged = []
        for s, e in hs:
            if merged and s - merged[-1][1] <= 50: merged[-1] = (merged[-1][0], max(merged[-1][1], e))
            else: merged.append((s, e))
        for s, e in merged:
            s2, e2 = max(0, s - pad), min(len(seq), e + pad)
            o.write(f">{c}:{s2}-{e2}\n{seq[s2:e2]}\n"); b.write(f"{c}\t{s2}\t{e2}\n"); n += 1; tot += e2 - s2
print(f"{n} candidate regions, {tot/1e6:.1f} Mbp")
PY
t3=$(date +%s)
if [[ $MODE == regions ]]; then
  $SS -g -3 -T $T -Q -n -z 11 -E 2 -Z 20000 -m 8C "$Q" "$O/regions.fa" 3 > "$O/pass2a.raw" 2>/dev/null
  t4=$(date +%s); echo "pass2a ssearch36 on regions: $((t4-t3)) s (regions built in $((t3-t2)) s), $(grep -vc '^#' "$O/pass2a.raw") HSPs"
  python3 "$B/sear_filter.py" "$O/pass2a.raw" --L $L 2>/dev/null | awk '{print $0"\tpass2a"}' > "$O/loci.bed"
  MASK="$O/regions.bed"
else
  # sear's own structure: filter pass-1 hits (fasta36 alignments), keep them, mask only those loci
  python3 "$B/sear_filter.py" "$O/pass1.raw" --L $L --rule $RULE --gap $GAP 2>/dev/null | awk '{print $0"\tpass1"}' > "$O/loci.bed"
  cut -f1-3 "$O/loci.bed" > "$O/plus.bed"; MASK="$O/plus.bed"
  t4=$(date +%s); echo "pass1 filtered ($RULE rule): $(wc -l < "$O/loci.bed") loci"
fi
if [[ $RES == 1 ]]; then
  python3 "$B/mask_hits.py" "$G" "$MASK" "$O/masked.fa" --pad 0 2>&1 | tail -1
  python3 - "$O/masked.fa" "$O/masked.win.fa" "$RES_WIN" "$RES_OVL" <<'PY'
import sys; sys.path.insert(0,'/home/user/knowlege/projects/sine-search-bench')
from bench import chunk_fasta, split_parts
chunk_fasta(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])); split_parts(sys.argv[2], 20000)
PY
  t5=$(date +%s)
  for p in "$O"/masked.win.fa.part*.fa; do $SS -g -3 -T $T -Q -n -z $RES_Z -E 2 -m 8C "$Q" "$p" 3; done > "$O/pass2b.raw" 2>/dev/null
  t6=$(date +%s); echo "pass2b ssearch36 on masked genome: $((t6-t5)) s (masking+windows $((t5-t4)) s), $(grep -vc '^#' "$O/pass2b.raw") HSPs"
  python3 "$B/sear_filter.py" "$O/pass2b.raw" --L $L --rule $RULE --gap $GAP 2>/dev/null | awk '{print $0"\tpass2b"}' >> "$O/loci.bed"
fi
echo "total $(( $(date +%s) - t0 )) s (incl. windowing $((t1-t0)) s); loci: $(wc -l < "$O/loci.bed") ($(grep -c pass2b "$O/loci.bed") from the residual pass)"
