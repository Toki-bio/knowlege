#!/usr/bin/env python3
"""Run several aligners on a simulated genome and score recall per divergence bin.

Usage: bench.py --dir SIMDIR [--tools ssearch36,blastn-w7,...] [--threads N]
Produces SIMDIR/results/<tool>.bed (chrom start end strand score) and
SIMDIR/results/summary.tsv.

Each tool entry: name -> (index_cmd or None, search_cmd, parser). Commands are
shell strings with {genome} {query} {chunks} {threads} {out} placeholders.
"""
import argparse, os, re, shlex, subprocess, sys, time
from collections import defaultdict

BIN = {
    "ssearch36": "/home/user/src/fasta36/bin/ssearch36",
    "fasta36": "/home/user/src/fasta36/bin/fasta36",
    "lalign36": "/home/user/src/fasta36/bin/lalign36",
    "blat": "/home/user/bin/blat",
    "bwa": "/home/user/src/bwa/bwa",
    "minimap2": "/home/user/src/minimap2/minimap2",
}


# ---------- parsers: return list of (chrom, start0, end, strand, score) ----------

def parse_blast_tab(path, chunked=False):
    hits = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        sid, ss, se, bits = f[1], int(f[8]), int(f[9]), float(f[11])
        strand = "+" if ss <= se else "-"
        s, e = min(ss, se) - 1, max(ss, se)
        if chunked:
            chrom, off = sid.rsplit(":", 1)[0], int(sid.rsplit(":", 1)[1].split("-")[0])
            s, e = s + off, e + off
            sid = chrom
        hits.append((sid, s, e, strand, bits))
    return hits


def parse_sam(path):
    hits = []
    for line in open(path):
        if line.startswith("@"):
            continue
        f = line.split("\t")
        flag, chrom, pos, cigar = int(f[1]), f[2], int(f[3]), f[5]
        if flag & 4 or chrom == "*":
            continue
        reflen = sum(int(n) for n, op in re.findall(r"(\d+)([MDN=X])", cigar))
        score = 0
        for t in f[11:]:
            if t.startswith("AS:i:"):
                score = int(t[5:])
        hits.append((chrom, pos - 1, pos - 1 + reflen, "-" if flag & 16 else "+", score))
    return hits


def parse_paf(path):
    hits = []
    for line in open(path):
        f = line.split("\t")
        if len(f) < 12:
            continue
        score = 0
        for t in f[12:]:
            if t.startswith("AS:i:"):
                score = int(t[5:])
        hits.append((f[5], int(f[7]), int(f[8]), f[4], score or int(f[9])))
    return hits


def parse_psl(path, merge_gap=50):
    """BLAT chains distant tiles as 'introns' (maxIntron 750 kb), so one PSL
    line can span several planted copies. Split into block runs instead."""
    hits = []
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 21 or not f[0].isdigit():
            continue
        sizes = [int(x) for x in f[18].rstrip(",").split(",")]
        tstarts = [int(x) for x in f[20].rstrip(",").split(",")]
        strand = f[8][-1] if len(f[8]) > 1 else f[8]
        run_s, run_e = tstarts[0], tstarts[0] + sizes[0]
        for st, sz in zip(tstarts[1:], sizes[1:]):
            if st - run_e > merge_gap:
                hits.append((f[13], run_s, run_e, strand, run_e - run_s))
                run_s = st
            run_e = st + sz
        hits.append((f[13], run_s, run_e, strand, run_e - run_s))
    return hits


def parse_nhmmer_tbl(path):
    hits = []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        chrom, af, at, strand, score = f[0], int(f[6]), int(f[7]), f[11], float(f[13])
        hits.append((chrom, min(af, at) - 1, max(af, at), strand, score))
    return hits


def parse_maf(path):
    """LAST MAF: 'a score=..' then two 's' lines; first s = genome (lastdb target)."""
    hits = []
    score = 0
    s_lines = []
    for line in open(path):
        if line.startswith("a "):
            m = re.search(r"score=(\d+)", line)
            score = int(m.group(1)) if m else 0
            s_lines = []
        elif line.startswith("s "):
            s_lines.append(line.split())
            if len(s_lines) == 2:
                name, start, aln_size, strand, src_size = s_lines[0][1], int(s_lines[0][2]), int(s_lines[0][3]), s_lines[0][4], int(s_lines[0][5])
                qstrand = s_lines[1][4]
                if strand == "-":  # LAST reports target on + normally; handle anyway
                    start = src_size - start - aln_size
                hits.append((name, start, start + aln_size, qstrand, score))
    return hits


# ---------- tool registry ----------
# (index_cmd, search_cmd, parser, needs_chunks)
T = {}
T["ssearch36"] = (None, f"{BIN['ssearch36']} -T {{threads}} -n -E 1e-3 -m 8C -z 3 {{query}} {{chunks}} > {{out}}", lambda p: parse_blast_tab(p, True))
T["ssearch36-fwd"] = (None, f"{BIN['ssearch36']} -T {{threads}} -n -E 1e-3 -m 8C -3 -z 3 {{query}} {{chunks}} > {{out}}", lambda p: parse_blast_tab(p, True))
T["ssearch36-parts"] = (None, "for p in {parts}; do " + BIN['ssearch36'] + " -T {threads} -n -E 1e-3 -m 8C -z 3 {query} $p; done > {out}", lambda p: parse_blast_tab(p, True))
T["fasta36-k4-parts"] = (None, "for p in {parts}; do " + BIN['fasta36'] + " -T {threads} -n -E 1e-3 -m 8C -z 3 {query} $p 4; done > {out}", lambda p: parse_blast_tab(p, True))
T["fasta36-k2"] = (None, f"{BIN['fasta36']} -T {{threads}} -n -E 1e-3 -m 8C -z 3 {{query}} {{chunks}} 2 > {{out}}", lambda p: parse_blast_tab(p, True))
T["fasta36-k4"] = (None, f"{BIN['fasta36']} -T {{threads}} -n -E 1e-3 -m 8C -z 3 {{query}} {{chunks}} 4 > {{out}}", lambda p: parse_blast_tab(p, True))
T["megablast"] = ("makeblastdb -in {genome} -dbtype nucl -out {genome}.db > /dev/null",
                  "blastn -task megablast -num_threads {threads} -query {query} -db {genome}.db -evalue 1e-3 -outfmt 6 -max_target_seqs 100000 -max_hsps 1000000 > {out}", parse_blast_tab)
T["dc-megablast"] = (None, "blastn -task dc-megablast -num_threads {threads} -query {query} -db {genome}.db -evalue 1e-3 -outfmt 6 -max_target_seqs 100000 -max_hsps 1000000 > {out}", parse_blast_tab)
T["blastn-w11"] = (None, "blastn -task blastn -num_threads {threads} -query {query} -db {genome}.db -evalue 1e-3 -outfmt 6 -max_target_seqs 100000 -max_hsps 1000000 > {out}", parse_blast_tab)
T["blastn-w7"] = (None, "blastn -task blastn -word_size 7 -num_threads {threads} -query {query} -db {genome}.db -evalue 1e-3 -outfmt 6 -max_target_seqs 100000 -max_hsps 1000000 > {out}", parse_blast_tab)
T["blastn-w7-nodust"] = (None, "blastn -task blastn -word_size 7 -dust no -soft_masking false -num_threads {threads} -query {query} -db {genome}.db -evalue 1e-3 -outfmt 6 -max_target_seqs 100000 -max_hsps 1000000 > {out}", parse_blast_tab)
T["blat-default"] = (None, f"{BIN['blat']} {{genome}} {{query}} -noHead {{out}} > /dev/null", parse_psl)
T["blat-sens"] = (None, f"{BIN['blat']} {{genome}} {{query}} -tileSize=8 -stepSize=4 -minScore=20 -minIdentity=50 -repMatch=1000000 -noHead {{out}} > /dev/null", parse_psl)
T["bwa-mem"] = (f"{BIN['bwa']} index {{genome}} 2> /dev/null", f"{BIN['bwa']} mem -t {{threads}} -a -T 20 {{genome}} {{query}} 2> /dev/null > {{out}}", parse_sam)
T["bwa-mem-k11"] = (None, f"{BIN['bwa']} mem -t {{threads}} -a -k 11 -c 100000 -T 20 -r 0.5 -y 100000 {{genome}} {{query}} 2> /dev/null > {{out}}", parse_sam)
T["minimap2-sr"] = (None, f"{BIN['minimap2']} -t {{threads}} -x sr -N 100000 --secondary=yes -c {{genome}} {{query}} 2> /dev/null > {{out}}", parse_paf)
T["minimap2-sens"] = (None, f"{BIN['minimap2']} -t {{threads}} -k 11 -w 3 -f 0 -n 1 -m 20 -s 20 -N 100000 -p 0.1 --secondary=yes -c {{genome}} {{query}} 2> /dev/null > {{out}}", parse_paf)
T["nhmmer"] = (None, "nhmmer --cpu {threads} -E 1e-3 --tblout {out} {query} {genome} > /dev/null", parse_nhmmer_tbl)
T["nhmmer-hmm"] = (None, "nhmmer --cpu {threads} -E 1e-3 --tblout {out} {query} {genome} > /dev/null", parse_nhmmer_tbl)
T["nhmmer-max"] = (None, "nhmmer --cpu {threads} --max -E 1e-3 --tblout {out} {query} {genome} > /dev/null", parse_nhmmer_tbl)
T["last-default"] = ("lastdb -P{threads} {genome}.lastdb {genome}", "lastal -P{threads} -E1e-3 {genome}.lastdb {query} > {out}", parse_maf)
T["last-near"] = ("lastdb -P{threads} -uNEAR {genome}.near {genome}", "lastal -P{threads} -E1e-3 {genome}.near {query} > {out}", parse_maf)
T["last-mam8"] = ("lastdb -P{threads} -uMAM8 {genome}.mam8 {genome}", "lastal -P{threads} -E1e-3 -m 100 -D 1e9 {genome}.mam8 {query} > {out}", parse_maf)
T["last-mam4"] = ("lastdb -P{threads} -uMAM4 {genome}.mam4 {genome}", "lastal -P{threads} -E1e-3 -m 100 -D 1e9 {genome}.mam4 {query} > {out}", parse_maf)
T["last-sens-m1000"] = (None, "lastal -P{threads} -E1e-3 -m 1000 -D 1e9 -r 2 -q 3 -a 5 -b 2 {genome}.lastdb {query} > {out}", parse_maf)
T["last-sens"] = (None, "lastal -P{threads} -E1e-3 -m 1000 -D 1e9 {genome}.lastdb {query} > {out}", parse_maf)


def split_parts(chunks, nrec=20000):
    """Split the chunk FASTA into files of <= nrec records (like sear's 20 parts);
    returns the list of part paths. ssearch36 keeps at most MAX_BEST=60,000 scan
    hits per run (weakest quarter dropped when full), so a genome-sized library
    must be searched in parts or weak copies are silently lost."""
    parts = []
    base = chunks + ".part"
    if os.path.exists(base + "0000.fa"):
        i = 0
        while os.path.exists(f"{base}{i:04d}.fa"):
            parts.append(f"{base}{i:04d}.fa"); i += 1
        return parts
    n, i, fh = 0, 0, None
    for line in open(chunks):
        if line.startswith(">"):
            if n % nrec == 0:
                if fh: fh.close()
                fh = open(f"{base}{i:04d}.fa", "w"); parts.append(f"{base}{i:04d}.fa"); i += 1
            n += 1
        fh.write(line)
    fh.close()
    return parts


def chunk_fasta(genome, out, size=5000, overlap=500):
    if os.path.exists(out):
        return
    from simulate import read_fasta
    with open(out, "w") as f:
        for name, seq in read_fasta(genome):
            for s in range(0, len(seq), size - overlap):
                e = min(len(seq), s + size)
                f.write(f">{name}:{s}-{e}\n{seq[s:e]}\n")
                if e == len(seq):
                    break


def load_truth(path):
    rows = []
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split("\t")
        rows.append((f[0], int(f[1]), int(f[2]), f[3], float(f[4]), f[5]))
    return rows


def score(hits, truth, min_frac=0.5):
    by_chrom = defaultdict(list)
    for h in hits:
        by_chrom[h[0]].append(h)
    for c in by_chrom:
        by_chrom[c].sort(key=lambda h: h[1])
    found = {}
    used = set()
    for chrom, s, e, name, div, strand in truth:
        best = 0.0
        for i, h in enumerate(by_chrom.get(chrom, [])):
            if h[2] <= s:
                continue
            if h[1] >= e:
                break
            ov = min(e, h[2]) - max(s, h[1])
            best = max(best, ov / (e - s))
            used.add((chrom, i))
        found[name] = best >= min_frac
    fp = sum(1 for c in by_chrom for i in range(len(by_chrom[c])) if (c, i) not in used)
    bins = defaultdict(lambda: [0, 0])
    for chrom, s, e, name, div, strand in truth:
        b = int(div * 100) // 5 * 5
        bins[b][1] += 1
        bins[b][0] += found[name]
    return found, fp, dict(sorted(bins.items()))


def run(cmd, cwd):
    t0 = time.time()
    r = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  ! exit {r.returncode}: {cmd}\n{r.stderr[-800:]}", file=sys.stderr)
    return time.time() - t0, r.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True)
    ap.add_argument("--tools", default=",".join(T))
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--query", default=None, help="override query fasta (e.g. an HMM for nhmmer)")
    ap.add_argument("--tag", default="", help="suffix for result names")
    a = ap.parse_args()
    d = os.path.abspath(a.dir)
    res = os.path.join(d, "results"); os.makedirs(res, exist_ok=True)
    genome, query = os.path.join(d, "genome.fa"), a.query or os.path.join(d, "query.fa")
    chunks = os.path.join(d, "genome.chunks.fa")
    chunk_fasta(genome, chunks)
    parts = " ".join(split_parts(chunks)) if any("-parts" in t for t in a.tools.split(",")) else ""
    truth = load_truth(os.path.join(d, "truth.bed"))
    summ_path = os.path.join(res, "summary.tsv")
    new = not os.path.exists(summ_path)
    with open(summ_path, "a") as summ:
        if new:
            summ.write("tool\tindex_s\tsearch_s\thits\tfp\trecall_all\t" + "\t".join(f"r{b}-{b+5}" for b in range(0, 40, 5)) + "\n")
        for name in a.tools.split(","):
            idx, cmd, parser = T[name]
            out = os.path.join(res, name + a.tag + ".raw")
            fmt = dict(genome=genome, query=query, chunks=chunks, parts=parts, threads=a.threads, out=out)
            ti = 0.0
            if idx:
                ti, rc = run(idx.format(**fmt), d)
            ts, rc = run(cmd.format(**fmt), d)
            hits = parser(out) if os.path.exists(out) else []
            with open(os.path.join(res, name + a.tag + ".bed"), "w") as f:
                for h in hits:
                    f.write("\t".join(map(str, h)) + "\n")
            found, fp, bins = score(hits, truth)
            rec = sum(found.values()) / len(found) if found else float("nan")
            row = [name + a.tag, f"{ti:.1f}", f"{ts:.1f}", str(len(hits)), str(fp), f"{rec:.3f}"] + [f"{bins[b][0]/bins[b][1]:.2f}" if b in bins and bins[b][1] else "-" for b in range(0, 40, 5)]
            summ.write("\t".join(row) + "\n"); summ.flush()
            print("\t".join(row))


if __name__ == "__main__":
    main()
