#!/usr/bin/env python3
"""Write the complement of BED intervals as FASTA records named chrom:start-end
(sear's minus bank: real sequence between found loci, no N masking).
Usage: complement_fasta.py genome.fa loci.bed out.fa [--min 100]"""
import sys, argparse
from collections import defaultdict
sys.path.insert(0, '/home/user/knowlege/projects/sine-search-bench')
from simulate import read_fasta
ap = argparse.ArgumentParser(); ap.add_argument("genome"); ap.add_argument("bed"); ap.add_argument("out"); ap.add_argument("--min", type=int, default=100)
a = ap.parse_args()
iv = defaultdict(list)
for l in open(a.bed):
    f = l.split("\t")
    if len(f) >= 3 and not l.startswith("#"): iv[f[0]].append((int(f[1]), int(f[2])))
n = tot = 0
with open(a.out, "w") as o:
    for name, seq in read_fasta(a.genome):
        pos = 0
        for s, e in sorted(iv.get(name, [])):
            if s - pos >= a.min:
                o.write(f">{name}:{pos}-{s}\n{seq[pos:s]}\n"); n += 1; tot += s - pos
            pos = max(pos, e)
        if len(seq) - pos >= a.min:
            o.write(f">{name}:{pos}-{len(seq)}\n{seq[pos:]}\n"); n += 1; tot += len(seq) - pos
print(f"{n} complement pieces, {tot/1e6:.1f} Mbp", file=sys.stderr)
