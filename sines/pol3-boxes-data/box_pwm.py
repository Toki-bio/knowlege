#!/usr/bin/env python3
"""box_pwm.py - Pol III A-box / B-box retention in SINE consensi vs real tRNA genes vs a shuffled null.

Standard library only (no numpy). Steps:
  1. Seed A/B box sites in every tRNA gene (loose regex inside a positional window; fallback = 11-mer with the
     fewest mismatches to the Galli consensus), build log-odds PWMs (pseudocount 0.5, uniform background),
     refine once (re-locate best PWM site in the window, rebuild). Leave-one-out self-scores give the tRNA score
     distribution; A->B start-to-start spacing is recorded (introns from tRNAscan-SE output if given).
  2. Score every SINE consensus (sense strand): best A anywhere; best B anywhere; best B downstream of the best A
     (spacing 20..120); Galli mismatches at the chosen 11-mers; a second A/B pair downstream (dimer check).
  3. Null: shuffled copies (esl-shuffle -d -N 20) scored the same way; per-sequence 95th percentile of the 20
     best-A (best-B) scores; also a pooled null over all shuffles.
  4. Rfam cmscan --tblout labels each consensus head: tRNA / 5S / 7SL(SRP) / other / none.

Usage:
  python3 box_pwm.py --trna hg19-tRNAs.fa --sines sines.fa --shuf sines.shuf.fa --cm cm.tbl
                     [--trnascan trnascan.out] [--controls Alu,B1,B2,ID] --out outdir
"""
import argparse
import collections
import math
import os
import re
import sys

GALLI_A = "TGGCNNAGTGG"   # Galli et al. 1981 A box (tRNA positions 8-18/19)
GALLI_B = "GGTTCGANNCC"   # B box (tRNA positions 52-62)
RE_A = re.compile(r"(?=(T[AG]G[CT]..A[AG].GG))")      # TRGYnnARnGG
RE_B = re.compile(r"(?=(GGTTC[AG]A..C[CT]))")         # GGTTCRAnnCY
L = 11
BASES = "ACGT"
IDX = {b: i for i, b in enumerate(BASES)}
PSEUDO = 0.5
BG = 0.25
SPACING_MIN, SPACING_MAX = 20, 120     # SINE A-start -> B-start window (tRNA ~44, Alu ~66)
HEAD_WIN = 40                          # a "head" A box must start within the first 40 nt


def read_fasta(path):
    recs, name = [], None
    with open(path) as fh:
        for l in fh:
            l = l.strip()
            if not l:
                continue
            if l.startswith(">"):
                name = l[1:].split()[0]
                recs.append([name, []])
            elif name is not None:
                recs[-1][1].append(l)
    return [(n, "".join(p).upper().replace("U", "T")) for n, p in recs]


def mismatches(kmer, cons):
    return sum(1 for a, b in zip(kmer, cons) if b != "N" and a != b)


def percentile(vals, q):
    """numpy-default (linear interpolation) percentile, q in 0..100."""
    v = sorted(vals)
    if not v:
        return float("nan")
    if len(v) == 1:
        return v[0]
    pos = (len(v) - 1) * q / 100.0
    lo, hi = int(math.floor(pos)), int(math.ceil(pos))
    return v[lo] + (v[hi] - v[lo]) * (pos - lo)


def median(vals):
    return percentile(vals, 50)


# ---------------------------------------------------------------- PWM
class PWM:
    def __init__(self, sites):
        self.n = len(sites)
        self.counts = [[0.0] * 4 for _ in range(L)]
        for s in sites:
            for i, c in enumerate(s):
                if c in IDX:
                    self.counts[i][IDX[c]] += 1
        self.rebuild()

    def rebuild(self):
        self.lo = []
        self.ic = []
        for i in range(L):
            tot = sum(self.counts[i]) + 4 * PSEUDO
            fr = [(self.counts[i][j] + PSEUDO) / tot for j in range(4)]
            self.lo.append([math.log2(f / BG) for f in fr])
            self.ic.append(sum(f * math.log2(f / BG) for f in fr))
        self.max_score = sum(max(r) for r in self.lo)
        self.consensus = "".join(BASES[r.index(max(r))] for r in self.lo)

    def score_kmer(self, kmer):
        s = 0.0
        for i, c in enumerate(kmer):
            s += self.lo[i][IDX[c]] if c in IDX else -2.0   # N penalised
        return s

    def scan(self, seq):
        """list of (score, start0) for every window."""
        return [(self.score_kmer(seq[i:i + L]), i) for i in range(len(seq) - L + 1)]

    def best(self, seq, lo=0, hi=None):
        """best (score, start0) with start0 in [lo, hi)."""
        if hi is None or hi > len(seq) - L + 1:
            hi = len(seq) - L + 1
        best = (float("-inf"), -1)
        for i in range(max(0, lo), hi):
            s = self.score_kmer(seq[i:i + L])
            if s > best[0]:
                best = (s, i)
        return best

    def loo_score(self, kmer):
        """score of a site that is part of the training set, with its own counts removed."""
        s = 0.0
        for i, c in enumerate(kmer):
            if c not in IDX:
                s += -2.0
                continue
            cnt = list(self.counts[i])
            cnt[IDX[c]] -= 1
            tot = sum(cnt) + 4 * PSEUDO
            s += math.log2(((cnt[IDX[c]] + PSEUDO) / tot) / BG)
        return s

    def table(self):
        out = ["pos\t" + "\t".join("lo_" + b for b in BASES) + "\t" + "\t".join("n_" + b for b in BASES) + "\tIC_bits"]
        for i in range(L):
            out.append("%d\t%s\t%s\t%.3f" % (i + 1, "\t".join("%.3f" % x for x in self.lo[i]),
                                            "\t".join("%d" % x for x in self.counts[i]), self.ic[i]))
        return "\n".join(out)


# ---------------------------------------------------------------- tRNA box seeding
def window_A(n):
    # A box starts at tRNA position 8 (1-based) in a canonical gene; allow 4..30 (1-based start)
    return 3, min(30, n - L + 1)


def window_B(n):
    # B box starts at position 52 in a 73-nt gene = 21 nt before the gene end (no CCA in GtRNAdb genomic seqs)
    return max(0, n - 34), n - L - 7 + 1


def seed_site(seq, regex, galli, lo, hi, anchor):
    """pick a regex match inside the window (fewest Galli mismatches, then closest to anchor);
    fallback: 11-mer in window with fewest Galli mismatches, ties closest to anchor.
    returns (start0, kmer, how)."""
    cands = []
    for m in regex.finditer(seq):
        i = m.start()
        if lo <= i < hi:
            cands.append((mismatches(m.group(1), galli), abs(i - anchor), i))
    if cands:
        cands.sort()
        i = cands[0][2]
        return i, seq[i:i + L], "regex"
    cands = [(mismatches(seq[i:i + L], galli), abs(i - anchor), i) for i in range(max(0, lo), hi)]
    if not cands:
        return -1, "", "none"
    cands.sort()
    i = cands[0][2]
    return i, seq[i:i + L], "mismatch"


def read_trnascan(path):
    """tRNAscan-SE tabular output -> {name: (intron_begin, intron_end)} (0,0 = no intron)."""
    d = {}
    if not path or not os.path.exists(path):
        return d
    for l in open(path):
        f = [x.strip() for x in l.rstrip("\n").split("\t")]
        if len(f) < 9 or not f[1].isdigit():
            continue
        try:
            ib, ie = int(f[6]), int(f[7])
        except ValueError:
            continue
        name = f[0]
        if name not in d or (ib and not d[name][0]):
            d[name] = (ib, ie)
    return d


def read_gtrnadb_detailed(fasta_path, detailed_path):
    """GtRNAdb hg19-tRNAs-detailed.out (tRNAscan-SE output keyed by chrN.trnaK) -> {gene name: (intron_begin, intron_end)}
    using the '(tRNAscan-SE ID: chr6.trna112)' tag of the GtRNAdb fasta headers."""
    d = {}
    if not detailed_path or not os.path.exists(detailed_path):
        return d
    introns = {}
    for l in open(detailed_path):
        f = [x.strip() for x in l.rstrip("\n").split("\t")]
        if len(f) < 9 or not f[1].isdigit():
            continue
        try:
            introns[f[0] + ".trna" + f[1]] = (int(f[6]), int(f[7]))
        except ValueError:
            continue
    for l in open(fasta_path):
        if l.startswith(">"):
            name = l[1:].split()[0]
            m = re.search(r"tRNAscan-SE ID:\s*(\S+?)\)", l)
            if m and m.group(1) in introns:
                d[name] = introns[m.group(1)]
    return d


def build_from_trnas(trnas, trnascan, outdir):
    rows = []
    sitesA, sitesB = [], []
    seeds = []
    for name, seq in trnas:
        n = len(seq)
        loA, hiA = window_A(n)
        loB, hiB = window_B(n)
        a0, ak, ahow = seed_site(seq, RE_A, GALLI_A, loA, hiA, 7)
        b0, bk, bhow = seed_site(seq, RE_B, GALLI_B, loB, hiB, n - 22)
        seeds.append((name, seq, a0, ak, ahow, b0, bk, bhow))
        if ak:
            sitesA.append(ak)
        if bk:
            sitesB.append(bk)
    pwmA0, pwmB0 = PWM(sitesA), PWM(sitesB)
    # one refinement pass: re-locate the best PWM site inside the window
    movedA = movedB = 0
    sitesA, sitesB, final = [], [], []
    for name, seq, a0, ak, ahow, b0, bk, bhow in seeds:
        n = len(seq)
        loA, hiA = window_A(n)
        loB, hiB = window_B(n)
        sa, a1 = pwmA0.best(seq, loA, hiA)
        sb, b1 = pwmB0.best(seq, loB, hiB)
        if a1 != a0:
            movedA += 1
        if b1 != b0:
            movedB += 1
        sitesA.append(seq[a1:a1 + L])
        sitesB.append(seq[b1:b1 + L])
        final.append((name, seq, a1, ahow, b1, bhow))
    pwmA, pwmB = PWM(sitesA), PWM(sitesB)
    # per-tRNA table with leave-one-out scores and spacing
    scoresA, scoresB, spac, spac_intron, spac_nointron = [], [], [], [], []
    with open(os.path.join(outdir, "trna_boxes.tsv"), "w") as fo:
        fo.write("tRNA\tlen\tintron\tA_start\tA_11mer\tA_seed\tA_mm_Galli\tA_score_LOO\t"
                 "B_start\tB_11mer\tB_seed\tB_mm_Galli\tB_score_LOO\tspacing_start_to_start\n")
        for name, seq, a1, ahow, b1, bhow in final:
            ak, bk = seq[a1:a1 + L], seq[b1:b1 + L]
            sa, sb = pwmA.loo_score(ak), pwmB.loo_score(bk)
            ib, ie = trnascan.get(name, (0, 0))
            intron = "%d-%d(%dnt)" % (ib, ie, abs(ie - ib) + 1) if ib else "no"
            if not trnascan:
                intron = "?"
            sp = b1 - a1
            scoresA.append(sa)
            scoresB.append(sb)
            spac.append(sp)
            (spac_intron if ib else spac_nointron).append(sp)
            fo.write("%s\t%d\t%s\t%d\t%s\t%s\t%d\t%.2f\t%d\t%s\t%s\t%d\t%.2f\t%d\n" % (
                name, len(seq), intron, a1 + 1, ak, ahow, mismatches(ak, GALLI_A), sa,
                b1 + 1, bk, bhow, mismatches(bk, GALLI_B), sb, sp))
    with open(os.path.join(outdir, "pwm_A.tsv"), "w") as fo:
        fo.write(pwmA.table() + "\n")
    with open(os.path.join(outdir, "pwm_B.tsv"), "w") as fo:
        fo.write(pwmB.table() + "\n")
    seedA = collections.Counter(x[3] for x in final)
    seedB = collections.Counter(x[5] for x in final)
    info = {
        "n_trna": len(final), "seedA": dict(seedA), "seedB": dict(seedB), "movedA": movedA, "movedB": movedB,
        "A_consensus": pwmA.consensus, "B_consensus": pwmB.consensus,
        "A_IC": sum(pwmA.ic), "B_IC": sum(pwmB.ic), "A_max": pwmA.max_score, "B_max": pwmB.max_score,
        "A_p5": percentile(scoresA, 5), "A_med": median(scoresA), "A_min": min(scoresA),
        "B_p5": percentile(scoresB, 5), "B_med": median(scoresB), "B_min": min(scoresB),
        "sp_med": median(spac), "sp_min": min(spac), "sp_max": max(spac), "sp_p5": percentile(spac, 5),
        "sp_p95": percentile(spac, 95),
        "sp_intron": (len(spac_intron), median(spac_intron) if spac_intron else float("nan"),
                      min(spac_intron) if spac_intron else 0, max(spac_intron) if spac_intron else 0),
        "sp_nointron": (len(spac_nointron), median(spac_nointron) if spac_nointron else float("nan"),
                        min(spac_nointron) if spac_nointron else 0, max(spac_nointron) if spac_nointron else 0),
        "spacing_hist": collections.Counter(spac),
    }
    return pwmA, pwmB, info


# ---------------------------------------------------------------- Rfam labels
def read_cmscan(path):
    """cmscan --tblout (default format): target name, accession, query name, ... ; returns {query: [(E, model, from, to, bits)]}"""
    d = collections.defaultdict(list)
    if not path or not os.path.exists(path):
        return d
    for l in open(path):
        if l.startswith("#"):
            continue
        f = l.split()
        if len(f) < 16:
            continue
        try:
            d[f[2]].append((float(f[15]), f[0], int(f[7]), int(f[8]), float(f[14]), f[9]))
        except ValueError:
            continue
    return d


def head_label(hits):
    if not hits:
        return "none", "-"
    hits = sorted(hits)
    E, model, a, b, bits, strand = hits[0]
    m = model.lower()
    if m == "trna" or m.startswith("trna"):
        lab = "tRNA"
    elif "5s_rrna" in m or m == "5s_rrna":
        lab = "5S"
    elif "srp" in m or "7sl" in m:
        lab = "7SL"
    else:
        lab = "other:" + model
    desc = "; ".join("%s %d-%d%s E=%.1e" % (h[1], h[2], h[3], "" if h[5] == "+" else "(-)", h[0]) for h in hits[:3])
    return lab, desc


# ---------------------------------------------------------------- SINE scoring
def score_seq(seq, pwmA, pwmB):
    """dict with best A, best B anywhere, best B downstream of A (spacing window), Galli mismatches."""
    r = {}
    if len(seq) < L:
        return None
    sa, a0 = pwmA.best(seq)
    sb_any, b_any = pwmB.best(seq)
    sb, b0 = pwmB.best(seq, a0 + SPACING_MIN, a0 + SPACING_MAX + 1)
    r.update(A_score=sa, A_pos=a0 + 1, A_11mer=seq[a0:a0 + L], A_mm_Galli=mismatches(seq[a0:a0 + L], GALLI_A),
             Bany_score=sb_any, Bany_pos=b_any + 1, Bany_11mer=seq[b_any:b_any + L],
             B_score=sb, B_pos=(b0 + 1) if b0 >= 0 else 0, B_11mer=seq[b0:b0 + L] if b0 >= 0 else "-",
             B_mm_Galli=mismatches(seq[b0:b0 + L], GALLI_B) if b0 >= 0 else -1,
             spacing=(b0 - a0) if b0 >= 0 else 0)
    # head region: best A starting within the first HEAD_WIN nt (SINE heads sit at the 5' end), B downstream of it
    sah, ah = pwmA.best(seq, 0, HEAD_WIN)
    sbh, bh = pwmB.best(seq, ah + SPACING_MIN, ah + SPACING_MAX + 1)
    r.update(Ahead_score=sah, Ahead_pos=ah + 1, Ahead_11mer=seq[ah:ah + L],
             Bhead_score=sbh, Bhead_pos=(bh + 1) if bh >= 0 else 0, Bhead_11mer=seq[bh:bh + L] if bh >= 0 else "-",
             spacing_head=(bh - ah) if bh >= 0 else 0)
    return r


def second_pair(seq, pwmA, pwmB, start0):
    """best A/B pair with A start >= start0 (dimer check); returns dict or None."""
    if start0 + 2 * L + SPACING_MIN > len(seq):
        return None
    sa, a0 = pwmA.best(seq, start0)
    if a0 < 0:
        return None
    sb, b0 = pwmB.best(seq, a0 + SPACING_MIN, a0 + SPACING_MAX + 1)
    if b0 < 0:
        return None
    return dict(A2_score=sa, A2_pos=a0 + 1, A2_11mer=seq[a0:a0 + L], B2_score=sb, B2_pos=b0 + 1,
                B2_11mer=seq[b0:b0 + L], spacing2=b0 - a0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trna", required=True)
    ap.add_argument("--sines", required=True)
    ap.add_argument("--shuf", required=True)
    ap.add_argument("--cm", default=None)
    ap.add_argument("--trnascan", default=None, help="tRNAscan-SE -o output run on the tRNA fasta (names = fasta names)")
    ap.add_argument("--trna_detailed", default=None, help="GtRNAdb hg19-tRNAs-detailed.out (intron coordinates by tRNAscan-SE ID)")
    ap.add_argument("--controls", default="Alu,AluY,B1,B2,ID,B1-dID,ID-Spe")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    trnas = read_fasta(a.trna)
    trnascan = read_trnascan(a.trnascan)
    if not trnascan and a.trna_detailed:
        trnascan = read_gtrnadb_detailed(a.trna, a.trna_detailed)
    pwmA, pwmB, info = build_from_trnas(trnas, trnascan, a.out)
    A_thr, B_thr = info["A_p5"], info["B_p5"]

    sines = read_fasta(a.sines)
    shuf = read_fasta(a.shuf)
    cm = read_cmscan(a.cm)

    # null per sequence
    nullA, nullB = collections.defaultdict(list), collections.defaultdict(list)
    for n, s in shuf:
        base = re.sub(r"-shuffled(-\d+)?$", "", n)
        r = score_seq(s, pwmA, pwmB)
        if r is None:
            continue
        nullA[base].append(r["A_score"])
        nullB[base].append(r["Bany_score"])
    pooledA = [x for v in nullA.values() for x in v]
    pooledB = [x for v in nullB.values() for x in v]
    pooled = dict(A_p95=percentile(pooledA, 95), A_p99=percentile(pooledA, 99), A_max=max(pooledA) if pooledA else float("nan"),
                  B_p95=percentile(pooledB, 95), B_p99=percentile(pooledB, 99), B_max=max(pooledB) if pooledB else float("nan"),
                  n=len(pooledA))

    # per-consensus table
    cols = ["name", "length", "head_Rfam", "Rfam_hits", "head_SINEBase", "A_score", "A_pos", "A_11mer", "A_mm_Galli",
            "B_score", "B_pos", "B_11mer", "B_mm_Galli", "spacing", "Bany_score", "Bany_pos", "Bany_11mer",
            "nullA_p95", "nullA_max", "nullB_p95", "nullB_max",
            "A_present_null", "A_present_tRNA", "B_present_null", "B_present_tRNA", "pair_present_tRNA",
            "A2_score", "A2_pos", "A2_11mer", "B2_score", "B2_pos", "B2_11mer", "spacing2", "dimer_tRNA",
            "Ahead_score", "Ahead_pos", "Ahead_11mer", "Bhead_score", "Bhead_pos", "Bhead_11mer", "spacing_head",
            "pair_head_tRNA"]
    rows = []
    for n, s in sines:
        r = score_seq(s, pwmA, pwmB)
        if r is None:
            continue
        lab, desc = head_label(cm.get(n, []))
        nA, nB = nullA.get(n, []), nullB.get(n, [])
        pA95 = percentile(nA, 95) if nA else float("nan")
        pB95 = percentile(nB, 95) if nB else float("nan")
        A_null = r["A_score"] > pA95 if nA else False
        B_null = r["Bany_score"] > pB95 if nB else False
        A_t = r["A_score"] >= A_thr
        B_t = r["Bany_score"] >= B_thr
        pair_t = A_t and r["B_pos"] > 0 and r["B_score"] >= B_thr
        pair_h = r["Ahead_score"] >= A_thr and r["Bhead_pos"] > 0 and r["Bhead_score"] >= B_thr
        row = dict(name=n, length=len(s), head_Rfam=lab, Rfam_hits=desc, head_SINEBase="", **r, pair_head_tRNA=int(pair_h),
                   nullA_p95=pA95, nullA_max=max(nA) if nA else float("nan"),
                   nullB_p95=pB95, nullB_max=max(nB) if nB else float("nan"),
                   A_present_null=int(A_null), A_present_tRNA=int(A_t), B_present_null=int(B_null),
                   B_present_tRNA=int(B_t), pair_present_tRNA=int(pair_t))
        # dimer check: second pair downstream of the first pair's B box (or of A if no B)
        first_end = (r["B_pos"] - 1 + L) if r["B_pos"] > 0 else (r["A_pos"] - 1 + L)
        p2 = second_pair(s, pwmA, pwmB, first_end) if pair_t else None
        if p2:
            row.update(p2)
            row["dimer_tRNA"] = int(p2["A2_score"] >= A_thr and p2["B2_score"] >= B_thr)
        else:
            row.update(A2_score=float("nan"), A2_pos=0, A2_11mer="-", B2_score=float("nan"), B2_pos=0, B2_11mer="-",
                       spacing2=0, dimer_tRNA=0)
        rows.append(row)

    def fmt(v):
        if isinstance(v, float):
            return "nan" if math.isnan(v) else "%.2f" % v
        return str(v)

    with open(os.path.join(a.out, "per_consensus.tsv"), "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for row in rows:
            fo.write("\t".join(fmt(row[c]) for c in cols) + "\n")

    # summary by head type
    def summarize(group, rows_):
        k = len(rows_)
        if not k:
            return None
        both_t = [r for r in rows_ if r["pair_present_tRNA"]]
        sp = [r["spacing"] for r in both_t]
        return [group, k,
                sum(r["A_present_tRNA"] for r in rows_), sum(r["B_present_tRNA"] for r in rows_), len(both_t),
                sum(r["A_present_null"] for r in rows_), sum(r["B_present_null"] for r in rows_),
                sum(1 for r in rows_ if r["A_present_null"] and r["B_present_null"]),
                sum(r["dimer_tRNA"] for r in rows_),
                median(sp) if sp else float("nan"), min(sp) if sp else 0, max(sp) if sp else 0,
                median([r["A_score"] for r in rows_]), median([r["Bany_score"] for r in rows_]),
                sum(r["pair_head_tRNA"] for r in rows_)]

    groups = collections.OrderedDict()
    for r in rows:
        g = r["head_Rfam"].split(":")[0]
        groups.setdefault(g, []).append(r)
    with open(os.path.join(a.out, "summary_by_head.tsv"), "w") as fo:
        fo.write("head_Rfam\tn\tA_tRNAthr\tB_tRNAthr\tpair_tRNAthr\tA_null95\tB_null95\tboth_null95\tdimers\t"
                 "spacing_median\tspacing_min\tspacing_max\tA_score_median\tB_score_median\tpair_head_tRNAthr\t"
                 "A_tRNAthr_pct\tB_tRNAthr_pct\tpair_tRNAthr_pct\n")
        for g in list(groups) + ["ALL"]:
            rr = rows if g == "ALL" else groups[g]
            s = summarize(g, rr)
            fo.write("\t".join(fmt(x) for x in s) + "\t%.1f\t%.1f\t%.1f\n" % (100.0 * s[2] / s[1], 100.0 * s[3] / s[1], 100.0 * s[4] / s[1]))
        fo.write("tRNA_genes(hg19)\t%d\t-\t-\t-\t-\t-\t-\t-\t%.1f\t%d\t%d\t%.2f\t%.2f\t-\t-\t-\t-\n" % (
            info["n_trna"], info["sp_med"], info["sp_min"], info["sp_max"], info["A_med"], info["B_med"]))

    with open(os.path.join(a.out, "A_without_B.tsv"), "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for r in rows:
            if r["A_present_tRNA"] and not r["B_present_tRNA"]:
                fo.write("\t".join(fmt(r[c]) for c in cols) + "\n")
    with open(os.path.join(a.out, "B_without_A.tsv"), "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for r in rows:
            if r["B_present_tRNA"] and not r["A_present_tRNA"]:
                fo.write("\t".join(fmt(r[c]) for c in cols) + "\n")
    with open(os.path.join(a.out, "dimers.tsv"), "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for r in rows:
            if r["dimer_tRNA"]:
                fo.write("\t".join(fmt(r[c]) for c in cols) + "\n")

    # controls
    want = [x.strip().lower() for x in a.controls.split(",") if x.strip()]
    with open(os.path.join(a.out, "controls.tsv"), "w") as fo:
        fo.write("\t".join(cols) + "\n")
        for r in rows:
            nl = r["name"].lower()
            if nl in want or any(nl.startswith(w) and len(nl) <= len(w) + 3 for w in want):
                fo.write("\t".join(fmt(r[c]) for c in cols) + "\n")

    # tRNA summary + null summary
    with open(os.path.join(a.out, "trna_summary.txt"), "w") as fo:
        fo.write("tRNA genes: %d\n" % info["n_trna"])
        fo.write("A seed method: %s ; B seed method: %s\n" % (info["seedA"], info["seedB"]))
        fo.write("sites moved in refinement: A %d, B %d\n" % (info["movedA"], info["movedB"]))
        fo.write("PWM A consensus %s  IC %.1f bits  max score %.2f\n" % (info["A_consensus"], info["A_IC"], info["A_max"]))
        fo.write("PWM B consensus %s  IC %.1f bits  max score %.2f\n" % (info["B_consensus"], info["B_IC"], info["B_max"]))
        fo.write("tRNA A score (leave-one-out): min %.2f  p5 %.2f  median %.2f\n" % (info["A_min"], info["A_p5"], info["A_med"]))
        fo.write("tRNA B score (leave-one-out): min %.2f  p5 %.2f  median %.2f\n" % (info["B_min"], info["B_p5"], info["B_med"]))
        fo.write("tRNA A->B spacing (start-to-start): median %.1f  p5 %.1f  p95 %.1f  min %d  max %d\n" % (
            info["sp_med"], info["sp_p5"], info["sp_p95"], info["sp_min"], info["sp_max"]))
        fo.write("  with intron (n=%d): median %.1f  range %d-%d\n" % info["sp_intron"])
        fo.write("  without intron (n=%d): median %.1f  range %d-%d\n" % info["sp_nointron"])
        fo.write("spacing histogram: %s\n" % " ".join("%d:%d" % kv for kv in sorted(info["spacing_hist"].items())))
        fo.write("null (pooled, %d shuffles): A p95 %.2f p99 %.2f max %.2f ; B p95 %.2f p99 %.2f max %.2f\n" % (
            pooled["n"], pooled["A_p95"], pooled["A_p99"], pooled["A_max"], pooled["B_p95"], pooled["B_p99"], pooled["B_max"]))
        fo.write("SINE consensi scored: %d ; with shuffles: %d\n" % (len(rows), len(nullA)))
        tl = open(os.path.join(a.out, "trna_boxes.tsv")).read().splitlines()[1:]
        ta = [float(l.split("\t")[7]) for l in tl]
        tb = [float(l.split("\t")[12]) for l in tl]
        fo.write("tRNA genes whose own A site scores above the pooled null p95 (%.2f): %d/%d ; above p99 (%.2f): %d/%d\n" % (
            pooled["A_p95"], sum(1 for x in ta if x > pooled["A_p95"]), len(ta),
            pooled["A_p99"], sum(1 for x in ta if x > pooled["A_p99"]), len(ta)))
        fo.write("tRNA genes whose own B site scores above the pooled null p95 (%.2f): %d/%d ; above p99 (%.2f): %d/%d\n" % (
            pooled["B_p95"], sum(1 for x in tb if x > pooled["B_p95"]), len(tb),
            pooled["B_p99"], sum(1 for x in tb if x > pooled["B_p99"]), len(tb)))
        fo.write("flags: *_tRNA = score >= 5th percentile of leave-one-out tRNA-gene site scores (A %.2f, B %.2f); "
                 "*_null = score > 95th percentile of the 20 shuffles of that consensus; "
                 "pair = A plus B 20..120 nt downstream of it; head = A starting within the first %d nt\n" % (A_thr, B_thr, HEAD_WIN))
    print(open(os.path.join(a.out, "trna_summary.txt")).read())
    print(open(os.path.join(a.out, "summary_by_head.tsv")).read())


if __name__ == "__main__":
    main()
