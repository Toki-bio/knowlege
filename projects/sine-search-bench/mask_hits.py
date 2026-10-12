#!/usr/bin/env python3
"""Mask hit intervals (BED: chrom start end ...) to N in a FASTA, with padding.

Usage: mask_hits.py genome.fa hits.bed out.fa [--pad 20]
Used for the two-pass experiment: mask what the fast pass found, then time the
sensitive pass on the residual genome.
"""
import sys, argparse
from collections import defaultdict

ap = argparse.ArgumentParser()
ap.add_argument("genome"); ap.add_argument("bed"); ap.add_argument("out")
ap.add_argument("--pad", type=int, default=20)
a = ap.parse_args()

iv = defaultdict(list)
for l in open(a.bed):
    f = l.split("\t")
    if len(f) < 3 or l.startswith("#"):
        continue
    iv[f[0]].append((int(f[1]), int(f[2])))

masked = total = 0
with open(a.genome) as fin, open(a.out, "w") as fout:
    name, buf = None, []
    def flush():
        global masked, total
        if name is None:
            return
        seq = bytearray("".join(buf), "ascii")
        total += len(seq)
        for s, e in iv.get(name, []):
            s, e = max(0, s - a.pad), min(len(seq), e + a.pad)
            seq[s:e] = b"N" * (e - s)
            masked += e - s
        fout.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            fout.write(seq[i:i+80].decode() + "\n")
    for line in fin:
        if line.startswith(">"):
            flush(); name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    flush()
print(f"masked {masked:,} of {total:,} bp ({100*masked/total:.2f}%) in {sum(len(v) for v in iv.values()):,} intervals", file=sys.stderr)
