#!/usr/bin/env python3
"""Score a BED of hits (chrom start end ...) against RepeatMasker L1 annotations.
Usage: score_l1.py hits.bed l1_all.bed [--minlen 4800] [--maxdiv 1000] [--chrom-map 21:chr21]
Truth = L1 rows with length >= minlen and milliDiv <= maxdiv. Found = >=50% of the element covered by hits.
Also reports hits overlapping no L1 row at all (candidate false positives)."""
import sys, argparse, bisect
from collections import defaultdict
ap = argparse.ArgumentParser(); ap.add_argument("hits"); ap.add_argument("truth"); ap.add_argument("--minlen", type=int, default=4800); ap.add_argument("--maxdiv", type=int, default=1000); ap.add_argument("--chrom-map", default="21:chr21"); ap.add_argument("--label", default="")
a = ap.parse_args()
cmap = dict(x.split(":") for x in a.chrom_map.split(",")) if a.chrom_map else {}
hits = defaultdict(list)
for l in open(a.hits):
    f = l.rstrip("\n").split("\t")
    if len(f) < 3 or l.startswith("#"): continue
    c = cmap.get(f[0], f[0]); hits[c].append((int(f[1]), int(f[2])))
for c in hits: hits[c].sort()
nh = sum(len(v) for v in hits.values())
truth = []; alll1 = defaultdict(list)
for l in open(a.truth):
    c, s, e, name, div, strand = l.rstrip("\n").split("\t")[:6]; s, e, div = int(s), int(e), int(div)
    alll1[c].append((s, e))
    if e - s >= a.minlen and div <= a.maxdiv: truth.append((c, s, e, name, div))
for c in alll1: alll1[c].sort()
def cover(c, s, e):
    L = hits.get(c, []); i = bisect.bisect_left(L, (s - 100000, 0)); cov = 0
    for hs, he in L[i:]:
        if hs >= e: break
        cov += max(0, min(e, he) - max(s, hs))
    return cov / (e - s)
found = [t for t in truth if cover(*t[:3]) >= 0.5]
# hits touching no L1 row
def touches(c, s, e):
    L = alll1.get(c, []); i = bisect.bisect_left(L, (s - 100000, 0))
    for ts, te in L[i:]:
        if ts >= e: break
        if te > s: return True
    return False
nol1 = sum(1 for c in hits for s, e in hits[c] if not touches(c, s, e))
print(f"{a.label}\thits={nh}\ttruth(L1>={a.minlen}bp,div<={a.maxdiv})={len(truth)}\tfound={len(found)}\trecall={len(found)/max(1,len(truth)):.3f}\thits_touching_no_L1={nol1}")
