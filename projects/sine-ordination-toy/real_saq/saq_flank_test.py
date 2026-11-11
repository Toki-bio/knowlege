"""Within-set unique-flank test for the 866 saq copies (method of Tal eri/LOG.md, 2026-08-21):
position-by-position identity of the flanks between every pair of copies (no alignment, which inflates identity);
a pair is 'related' if identity > 45 % (background: 0.1 % of random 70 bp genomic window pairs exceed 45 %).
Flanks as stored in the Tal rand100 rows: first 50 bp = left, last 70 bp = right.
Caveat from eri/LOG.md: the right 70 bp can still be SINE tail; there the test needed 100-150 bp more downstream.
Also reported: right flank without its first 30 bp (positions 31-70), as a partial guard against tail.
Output: saq_flank_pairs.tsv (copy, left/right partners > 45 %, max identity), summary to stdout + saq_flank_summary.json
"""
import json, re, numpy as np
seqs = {}; h = None
for l in open("saq_rand100_copies.fa"):
    l = l.strip()
    if l.startswith(">"): h = l[1:]; seqs[h] = ""
    else: seqs[h] += l.upper()
names = [n for n, s in seqs.items() if len(s) - (lambda a, b: b - a)(*map(int, re.search(r":(\d+)-(\d+)\(", n).groups())) == 120]
code = {"A": 0, "C": 1, "G": 2, "T": 3}
def arr(parts): return np.array([[code.get(c, 9) for c in p] for p in parts])
Lf = arr([seqs[n][:50] for n in names]); Rf = arr([seqs[n][-70:] for n in names]); Rf2 = Rf[:, 30:]
def pair_ident(F):
    n = len(F); I = np.zeros((n, n))
    valid = F < 4
    for i in range(n):
        both = valid[i] & valid
        I[i] = ((F[i] == F) & both).sum(1) / np.maximum(both.sum(1), 1)
    np.fill_diagonal(I, 0); return I
lab = np.array([n.split("|")[1] for n in names])
out = {}
res = {}
for side, F in (("left50", Lf), ("right70", Rf), ("right31_70", Rf2)):
    I = pair_ident(F); iu = np.triu_indices(len(F), 1)
    res[side] = I
    out[side] = dict(mean_pair_identity=float(I[iu].mean()), frac_pairs_gt45=float((I[iu] > 0.45).mean()),
                     copies_with_partner_gt45=int((I > 0.45).any(1).sum()),
                     per_subfamily_copies_with_partner={s: int(((I > 0.45).any(1) & (lab == s)).sum()) for s in sorted(set(lab))})
twin = (res["left50"] > 0.45) & (res["right70"] > 0.45)
out["copies_with_same_partner_both_sides_gt45"] = int(twin.any(1).sum())
with open("saq_flank_pairs.tsv", "w") as f:
    f.write("copy\tsubfamily\tleft_partners_gt45\tright_partners_gt45\tboth_sides_partners\tmax_left\tmax_right\n")
    for i, n in enumerate(names):
        f.write(f"{n}\t{lab[i]}\t{int((res['left50'][i] > .45).sum())}\t{int((res['right70'][i] > .45).sum())}\t{int(twin[i].sum())}\t"
                f"{res['left50'][i].max():.2f}\t{res['right70'][i].max():.2f}\n")
json.dump(out, open("saq_flank_summary.json", "w"), indent=1); print(json.dumps(out, indent=1))
