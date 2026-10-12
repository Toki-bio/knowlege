#!/usr/bin/env python3
"""Apply sear's post-filter to a FASTA36 -m 8C output over named windows.

Two variants of the length rule:
  hsp : max single-HSP alignment length >= slf*L      (sear at 85421e7, col 4 = max)
  span: merged-span length >= slf*L, HSPs merged within --gap bp
Identity rule: max %id over the merged HSPs >= hom. Prints loci as BED.
"""
import argparse, collections, sys
ap = argparse.ArgumentParser()
ap.add_argument("raw"); ap.add_argument("--L", type=int, default=357); ap.add_argument("--slf", type=float, default=0.8)
ap.add_argument("--hom", type=float, default=65); ap.add_argument("--rule", choices=["hsp", "span"], default="hsp"); ap.add_argument("--gap", type=int, default=0)
a = ap.parse_args()
hits = collections.defaultdict(list)
for l in open(a.raw):
    if l.startswith("#") or not l.strip(): continue
    f = l.rstrip("\n").split("\t"); chrom, rng = f[1].rsplit(":", 1); off = int(rng.split("-")[0])
    s, e = sorted((int(f[8]), int(f[9]))); hits[chrom].append((s - 1 + off, e + off, int(f[3]), float(f[2]), "+" if int(f[7]) >= int(f[6]) else "-"))
n = 0
for c in sorted(hits):
    hs = sorted(hits[c]); merged = []
    for s, e, ln, pid, st in hs:
        if merged and s - merged[-1][1] <= a.gap:
            ps, pe, pl, pp, pst = merged[-1]; merged[-1] = (ps, max(pe, e), max(pl, ln), max(pp, pid), pst)
        else: merged.append((s, e, ln, pid, st))
    for s, e, ln, pid, st in merged:
        length = ln if a.rule == "hsp" else e - s
        if length >= a.slf * a.L and pid >= a.hom:
            print(f"{c}\t{s}\t{e}\t{st}\t{pid}"); n += 1
print(f"{n} loci ({a.rule} rule, gap {a.gap})", file=sys.stderr)
