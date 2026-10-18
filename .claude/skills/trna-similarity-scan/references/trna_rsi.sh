#!/bin/bash
# tRNA / 5S / 7SL origin of the rsi consensuses and of the composites' parts (2026-10-18).
# Queries: the 21 consensuses of rsi_main_20261018 + every stretch of the 7 composites from composite_units.tsv.
# Searches: (1) ssearch36 vs 419 human tRNA genes (GtRNAdb hg19) + Rfam-free 5S/7SL from the partner library,
#           SINEBase rule >= 60 % over >= 60 nt; best and 2nd-best tRNA isotype; null = 20 dinucleotide shuffles per query
#           (esl-shuffle -d), best shuffled identity*length reported;
#           (2) Infernal cmscan vs Rfam (E <= 1), (3) tRNAscan-SE -E -H --detail (isotype scores), (4) A/B box scan.
set -uo pipefail
source ~/miniforge3/etc/profile.d/conda.sh
OUT=/home/toki/rhin/rsi_main_20261018/trna; mkdir -p $OUT; cd $OUT || exit 1
export TMPDIR=$HOME/tmp
R=$(cat /home/toki/rhin/rsi_main_20261018/RUN)
conda activate sinederella
python3 - "$R" <<'PY'
import sys, csv, glob
R = sys.argv[1]
seq, n = {}, None
for l in open(R + "/consensuses.clean.fa"):
    l = l.strip()
    if l.startswith(">"): n = l[1:].split()[0]; seq[n] = ""
    else: seq[n] += l.upper()
cu = sorted(glob.glob(R + "/step2/step2_output*/composite_units.tsv"))[-1]
with open("queries.fa", "w") as fo:
    for k, s in seq.items():
        fo.write(">%s\n%s\n" % (k, s))
    for r in csv.DictReader(open(cu), delimiter="\t"):
        if r["kind"] != "composite": continue
        a, b = int(r["start"]), int(r["end"])
        if b - a + 1 < 30: continue
        fo.write(">%s:%d-%d|%s\n%s\n" % (r["consensus"], a, b, r["unit"], seq[r["consensus"]][a - 1:b]))
PY
echo "queries: $(grep -c '>' queries.fa)"
# reference: human tRNA genes + 5S / 7SL entries of the partner library
grep -A1 -i -E "^>x\.(5S|7SL|SRP)" /home/toki/refs/partnerlib/partnerlib_9397.fa | grep -v -- "^--$" > rna_5s7sl.fa
cat /home/toki/refs/smallrna/hg19-tRNAs.fa rna_5s7sl.fa > ref.fa
echo "reference: $(grep -c '>' /home/toki/refs/smallrna/hg19-tRNAs.fa) tRNA + $(grep -c '>' rna_5s7sl.fa) 5S/7SL"
ssearch36 -3 -m 8 -E 10 -T 16 queries.fa ref.fa > ss.m8 2>/dev/null
# shuffled null: 20 dinucleotide shuffles of each query
esl-shuffle -d -N 20 --seed 42 queries.fa > shuf.fa 2>/dev/null || python3 - <<'PY'
import random
random.seed(42)
seqs = []; n = None
for l in open("queries.fa"):
    l = l.strip()
    if l.startswith(">"): n = l[1:]; seqs.append([n, ""])
    else: seqs[-1][1] += l
with open("shuf.fa", "w") as fo:
    for n, s in seqs:
        for k in range(20):
            x = list(s); random.shuffle(x)
            fo.write(">%s-shuffled-%d\n%s\n" % (n, k + 1, "".join(x)))
PY
ssearch36 -3 -m 8 -E 10 -T 16 shuf.fa ref.fa > shuf.m8 2>/dev/null
# Infernal vs Rfam, tRNAscan-SE (rnatools env)
conda activate rnatools
cmscan --cpu 16 -E 1 --tblout cm.tbl -o /dev/null /home/toki/refs/smallrna/Rfam.cm queries.fa
tRNAscan-SE -E -H --detail -X 10 -q -o trnascan.out -f trnascan.ss queries.fa > /dev/null 2>&1 || tRNAscan-SE -E -H -X 10 -q -o trnascan.out queries.fa > /dev/null 2>&1
echo "=== done $(date)"
