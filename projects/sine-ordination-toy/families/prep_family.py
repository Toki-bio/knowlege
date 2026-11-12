"""Prepare one family from the Tal repository for sine_inputs_v4 / family_inputs_v4.

Input : <tal_dir>/alignments/<prefix>_<sub>_rand100.aln.fa for each listed subfamily. Each file is an alignment of
        the subfamily consensus (row named <sub>) and up to 100 random copies with flanks.
Body  : the alignment columns from the first to the last base of the consensus row; flanks = what lies outside.
Flank gate (method of Tal eri/LOG.md, as in real_saq/saq_flank_test.py): position-by-position identity of the
        50 bp left and 70 bp right flanks between every two copies, no alignment; a copy is DROPPED if with some
        other copy it is > 60 % identical on both sides. Copies with a flank shorter than 30 bp are dropped too
        (their flanks cannot be tested).
Output: <out>_body.fa (header coord|sub), <out>_consensuses.fa, <out>_prep.json
Usage : python prep_family.py TAL_DIR PREFIX OUT sub1 sub2 ...
"""
import sys, json, re, numpy as np
tal, prefix, out, subs = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
def fa(p):
    d = {}; h = None
    for l in open(p, encoding="utf-8"):
        l = l.strip()
        if l.startswith(">"): h = l[1:]; d[h] = ""
        elif h is not None: d[h] += l
    return d
rows = []; cons = {}; info = dict(subfamilies=subs, per_subfamily={})
for s in subs:
    d = fa(f"{tal}/alignments/{prefix}_{s}_rand100.aln.fa")
    if s not in d: sys.exit(f"no consensus row {s}")
    c = d[s]; pos = [i for i, x in enumerate(c) if x != "-"]; a, b = pos[0], pos[-1] + 1
    cons[s] = c.replace("-", "").upper(); n = 0
    for h, q in d.items():
        if not re.search(r":\d+-\d+\(", h): continue
        body = q[a:b].replace("-", "").upper(); L = q[:a].replace("-", "").upper(); R = q[b:].replace("-", "").upper()
        rows.append((h.split()[0], s, body, L[-50:], R[:70])); n += 1
    info["per_subfamily"][s] = dict(copies=n, consensus_len=len(cons[s]))
code = {"A": 0, "C": 1, "G": 2, "T": 3}
ok = np.array([len(r[3]) >= 30 and len(r[4]) >= 30 and len(r[2]) >= 30 for r in rows])
def arr(parts, n, left):
    M = np.full((len(parts), n), 9)
    for i, p in enumerate(parts):
        v = [code.get(x, 9) for x in p]
        if left: M[i, n - len(v):] = v          # left flank: aligned at the element side
        else: M[i, :len(v)] = v
    return M
Lf = arr([r[3] for r in rows], 50, True); Rf = arr([r[4] for r in rows], 70, False)
def ident(F):
    n = len(F); I = np.zeros((n, n)); V = F < 4
    for i in range(n):
        both = V[i] & V; m = both.sum(1); I[i] = np.where(m >= 30, ((F[i] == F) & both).sum(1) / np.maximum(m, 1), 0)
    np.fill_diagonal(I, 0); return I
IL, IR = ident(Lf), ident(Rf); twin = (IL > 0.6) & (IR > 0.6); twin &= ok[:, None] & ok[None, :]
shared = twin.any(1); keep = ok & ~shared
info.update(copies_in=len(rows), dropped_short_flank_or_body=int((~ok).sum()), dropped_shared_flanks=int(shared.sum()), copies_kept=int(keep.sum()),
            mean_pair_identity_left=float(IL[np.ix_(ok, ok)].mean()), mean_pair_identity_right=float(IR[np.ix_(ok, ok)].mean()),
            kept_per_subfamily={s: int(sum(1 for r, k in zip(rows, keep) if k and r[1] == s)) for s in subs})
with open(out + "_body.fa", "w") as f:
    for r, k in zip(rows, keep):
        if k: f.write(f">{r[0]}|{r[1]}\n{r[2]}\n")
with open(out + "_consensuses.fa", "w") as f:
    for s in subs: f.write(f">{s}\n{cons[s]}\n")
json.dump(info, open(out + "_prep.json", "w"), indent=1); print(json.dumps(info, indent=1))
