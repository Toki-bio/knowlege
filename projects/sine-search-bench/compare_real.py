#!/usr/bin/env python3
"""Compare hit sets on a real genome where no ground truth exists.

Reference = one tool's BED (default ssearch36). For every other tool, report
how many reference hits it recovers (>=50% overlap), binned by the reference
hit's percent identity when available (ssearch36 -m 8C column 3), and how many
hits it has that the reference lacks ("extra").

Usage: compare_real.py RESULTS_DIR [--ref ssearch36] [--ref-raw ssearch36.raw]
"""
import argparse, os, sys
from collections import defaultdict

def load_bed(p, merge_gap=50):
    """Load hits and merge overlapping / near-adjacent ones on a chromosome so
    that a copy reported as two HSPs counts once (strand ignored)."""
    d = defaultdict(list)
    for l in open(p):
        f = l.split("\t")
        d[f[0]].append((int(f[1]), int(f[2]), f[3], float(f[4])))
    for c in d:
        d[c].sort()
        merged = []
        for s, e, strand, sc in d[c]:
            if merged and s - merged[-1][1] <= merge_gap:
                ps, pe, pstr, psc = merged[-1]
                merged[-1] = (ps, max(pe, e), pstr, max(psc, sc))
            else:
                merged.append((s, e, strand, sc))
        d[c] = merged
    return d

def load_ref_ident(raw, chunked=True):
    """ssearch36 -m 8C: identity per hit, keyed by (chrom,start,end)."""
    ident = {}
    for l in open(raw):
        if l.startswith("#") or not l.strip():
            continue
        f = l.rstrip("\n").split("\t")
        sid, ss, se = f[1], int(f[8]), int(f[9])
        s, e = min(ss, se) - 1, max(ss, se)
        if chunked:
            chrom, off = sid.rsplit(":", 1)[0], int(sid.rsplit(":", 1)[1].split("-")[0])
            s, e, sid = s + off, e + off, chrom
        ident[(sid, s, e)] = float(f[2])
    return ident

def overlap_fraction(ref, hits):
    """For each ref interval, max overlap fraction by any hit; and set of hit indices used."""
    import bisect
    out = {}
    used = set()
    for c, ivs in ref.items():
        hs = hits.get(c, [])
        starts = [h[0] for h in hs]
        for s, e, strand, sc in ivs:
            i = bisect.bisect_left(starts, s - 5000)
            best = 0.0
            while i < len(hs) and hs[i][0] < e:
                hs_, he_ = hs[i][0], hs[i][1]
                ov = min(e, he_) - max(s, hs_)
                if ov > 0:
                    best = max(best, ov / (e - s)); used.add((c, i))
                i += 1
            out[(c, s, e)] = best
    return out, used

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir"); ap.add_argument("--ref", default="ssearch36"); ap.add_argument("--min-frac", type=float, default=0.5)
    a = ap.parse_args()
    ref = load_bed(os.path.join(a.dir, a.ref + ".bed"))
    raw = os.path.join(a.dir, a.ref + ".raw")
    ident = load_ref_ident(raw) if a.ref.startswith(("ssearch", "fasta")) and os.path.exists(raw) else {}
    nref = sum(len(v) for v in ref.values())
    bins = [(95, 101), (90, 95), (85, 90), (80, 85), (75, 80), (70, 75), (60, 70), (0, 60)]
    refbin = defaultdict(int)
    for k, pid in ident.items():
        for lo, hi in bins:
            if lo <= pid < hi: refbin[(lo, hi)] += 1
    print(f"reference {a.ref}: {nref} hits; identity bins: " + ", ".join(f"{lo}-{hi}:{refbin[(lo,hi)]}" for lo, hi in bins))
    print("tool\thits\trecall_of_ref\textra_hits\t" + "\t".join(f"r{lo}-{hi}" for lo, hi in bins))
    for f in sorted(os.listdir(a.dir)):
        if not f.endswith(".bed") or f == a.ref + ".bed":
            continue
        hits = load_bed(os.path.join(a.dir, f))
        nh = sum(len(v) for v in hits.values())
        fr, used = overlap_fraction(ref, hits)
        rec = sum(1 for v in fr.values() if v >= a.min_frac) / max(1, nref)
        extra = nh - len(used)
        row = [f[:-4], str(nh), f"{rec:.3f}", str(extra)]
        for lo, hi in bins:
            keys = [k for k, pid in ident.items() if lo <= pid < hi and k in fr]
            row.append(f"{sum(1 for k in keys if fr[k] >= a.min_frac)/len(keys):.2f}" if keys else "-")
        print("\t".join(row))

if __name__ == "__main__":
    main()
