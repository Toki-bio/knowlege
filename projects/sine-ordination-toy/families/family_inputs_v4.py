"""One real family: K2P over all columns against positions selected by covariation (method: ../sine_inputs_v4.py).

Input : MAFFT alignment of the subfamily consensuses (rows without '|') and the copy bodies (rows 'coord|label').
Output: <P>_v4.json, <P>_v4_features.tsv, <P>_v4_crosstabs.txt, <P>_v4.png and two alignments for reading:
        <P>_groups_selected.aln.fa  consensuses, then copies ordered by nested group; ONLY the selected columns
        <P>_groups_full.aln.fa      same rows and order; all columns covered by at least half the copies
        Copy names: g<group, substitutions only>_i<group, substitutions + indels>|<published label>|<coordinates>
Usage : python family_inputs_v4.py ALN PREFIX
"""
import sys, json, numpy as np
sys.path.insert(0, ".."); sys.path.insert(0, "/home/toki/ordination_toy")
aln, P = sys.argv[1], sys.argv[2]; sys.argv = sys.argv[:1]
import sine_inputs_v4 as s4
import ordination_toy_v3 as v3
from sklearn.metrics import adjusted_rand_score
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

SEED = 1
names, seqs = [], []
for line in open(aln):
    line = line.strip()
    if line.startswith(">"): names.append(line[1:]); seqs.append([])
    elif line: seqs[-1].append(line.upper())
seqs = ["".join(s) for s in seqs]; code = {"A": 0, "C": 1, "G": 2, "T": 3}
A = np.array([[code.get(c, -1) for c in s] for s in seqs]); iscopy = np.array(["|" in n for n in names])
X = A[iscopy]; cname = [n for n, c in zip(names, iscopy) if c]; cseq = [s for s, c in zip(seqs, iscopy) if c]
Cn = [n for n, c in zip(names, iscopy) if not c]; Cseq = [s for s, c in zip(seqs, iscopy) if not c]; C = A[~iscopy].copy(); C[C < 0] = s4.GAPS
lab_s = np.array([n.split("|")[1] for n in cname]); subs = [c for c in Cn if c in set(lab_s)] or sorted(set(lab_s)); lab = np.array([subs.index(s) for s in lab_s])
N, L = X.shape; S = s4.code_states(X); LET = "ACGT-"
out = dict(copies=int(N), align_cols=int(L), subfamilies=subs, copies_per_subfamily={s: int((lab_s == s).sum()) for s in subs})
# consensus relatedness: identity over columns where both consensuses have a base
ci = {}
for a in range(len(Cn)):
    for b in range(a + 1, len(Cn)):
        m = (C[a] != s4.GAPS) & (C[b] != s4.GAPS); ci[f"{Cn[a]}~{Cn[b]}"] = [int(m.sum()), round(float((C[a][m] == C[b][m]).mean()), 3) if m.sum() else None]
out["consensus_pairs_shared_cols_identity"] = ci
xt = open(f"{P}_v4_crosstabs.txt", "w")
def crosstab(title, cl):
    xt.write(f"## {title}\ngroup\t" + "\t".join(subs) + "\ttotal\n")
    for c in sorted(set(cl)):
        row = [int(((cl == c) & (lab == i)).sum()) for i in range(len(subs))]; xt.write(f"{c}\t" + "\t".join(map(str, row)) + f"\t{sum(row)}\n")
    xt.write("\n")
pics = {}
def arm(name, D):
    sc, U, T, hd = s4.evaluate(D, lab, SEED, len(subs)); out[name] = sc; pics[name] = (U, hd); crosstab(f"{name}: HDBSCAN on UMAP x published", hd)
    print(name, json.dumps(sc), flush=True)

cover = (X >= 0).mean(0); keep = cover >= 0.5
cons = np.array([np.bincount(c[c >= 0], minlength=4).argmax() if (c >= 0).any() else -1 for c in X.T])
cg = np.zeros(L, bool); k = np.where(keep)[0]; ck = cons[k]; cgk = (ck[:-1] == 1) & (ck[1:] == 2); cg[k[:-1][cgk]] = True; cg[k[1:][cgk]] = True
out.update(cols_half_covered=int(keep.sum()), cpg_cols=int(cg.sum()))
Dk = v3.k2p_pd(X[:, keep & ~cg]); out["saturated_pairs_K2P"] = int((Dk >= v3.CLIP - 1e-9).sum() // 2); arm("K2P CpG-masked", Dk)

res = {}
for tag, indel in (("sub", False), ("sub+indel", True)):
    ft = s4.build_features(S, indel=indel); w, c = s4.cov_weights(ft, S, SEED); sel = w > 0; res[tag] = (ft, w, sel)
    out[f"selection {tag}"] = dict(features=int(ft["Fm"].shape[1]), selected=int(sel.sum()), selected_indel=int((sel & (ft["kind"] != "sub")).sum()),
                                   threshold_z2=c, selected_on_CpG_cols=int(cg[ft["col"][sel]].sum()))
    if tag == "sub": arm("all characters", s4.wham(ft["Fm"], np.ones(ft["Fm"].shape[1])))
    if sel.sum() >= 2: arm(f"selected {tag}", s4.wham(ft["Fm"][:, sel], np.ones(int(sel.sum()))))
ft, w, sel = res["sub+indel"]

cvar = np.array([len(set(C[:, j].tolist())) > 1 for j in range(L)])
def cons_hit(j, kind, state):
    cs = C[:, j]
    if not cvar[j]: return False
    return bool((cs != s4.GAPS).any() and (cs == s4.GAPS).any()) if kind == "ins" else bool((cs == state).any())
with open(f"{P}_v4_features.tsv", "w") as f:
    f.write("column\tkind\tstate\tfreq\tweight\tselected\tCpG_col\tconsensuses_differ_here\tstate_is_a_consensus_state\tconsensus_states_" + ",".join(Cn) + "\t" + "\t".join("freq_" + s for s in subs) + "\n")
    nhit = 0
    for i in np.argsort(-w):
        j = int(ft["col"][i]); hit = cons_hit(j, ft["kind"][i], ft["state"][i]); nhit += bool(sel[i] and hit)
        with np.errstate(all="ignore"): fr = [np.nanmean(ft["Fm"][lab == q, i]) for q in range(len(subs))]
        f.write(f"{j}\t{ft['kind'][i]}\t{LET[ft['state'][i]] if ft['state'][i] >= 0 else 'base'}\t{np.nanmean(ft['Fm'][:, i]):.3f}\t{w[i]:.1f}\t{int(sel[i])}\t{int(cg[j])}\t{int(cvar[j])}\t{int(hit)}\t{''.join(LET[x] for x in C[:, j])}\t" + "\t".join(f"{x:.2f}" for x in fr) + "\n")
fcols = set(int(j) for j in ft["col"]); selcols = set(int(j) for j in ft["col"][sel]); cvf = set(j for j in fcols if cvar[j])
out["consensus check (sub+indel)"] = dict(selected_characters=int(sel.sum()), selected_that_are_a_consensus_difference=int(nhit),
                                         cols_where_consensuses_differ_and_a_character_exists=len(cvf), of_those_selected=len(cvf & selcols),
                                         selected_cols=len(selcols), selected_cols_where_consensuses_do_not_differ=len(selcols - set(np.where(cvar)[0].tolist())))
cd = np.array([cvar[j] for j in ft["col"]])
if cd.sum() >= 2: arm("consensus-difference columns (NOT label-free)", s4.wham(ft["Fm"][:, cd], np.ones(int(cd.sum()))))

rng = np.random.default_rng(7); halves = []
for rep in range(5):
    p = rng.permutation(N); h = []
    for part in (p[:N // 2], p[N // 2:]):
        fh = s4.build_features(S[part], indel=True); wh, _ = s4.cov_weights(fh, S[part], SEED + rep)
        h.append(set((int(a), str(b), int(c_)) for a, b, c_ in zip(fh["col"][wh > 0], fh["kind"][wh > 0], fh["state"][wh > 0])))
    halves.append(dict(a=len(h[0]), b=len(h[1]), both=len(h[0] & h[1])))
out["split-half selection (5 random splits)"] = halves

nl = {}; dcols = set()
for tag, indel in (("sub", False), ("sub+indel", True)):
    log = []; lv = s4.nested(S, SEED, indel=indel, log=log); nl[tag] = s4.leaves_to_labels(lv, N)
    for n_ in log:
        for v in n_.get("diag_positions_between_parts", {}).values(): dcols |= set(v)
    out[f"nested {tag}"] = dict(groups=len(lv), sizes=[int(len(v)) for _, v in sorted(lv.items())], ari_vs_published=float(adjusted_rand_score(lab, nl[tag])), log=log)
    crosstab(f"nested {tag}: groups x published (group order = {sorted(lv)})", nl[tag])
xt.close()

# alignments for reading
order = np.lexsort((lab, nl["sub+indel"], nl["sub"])); showc = sorted(selcols | dcols | set(int(j) for j in res["sub"][0]["col"][res["sub"][2]]))
fullc = [j for j in range(L) if keep[j] or (C[:, j] != s4.GAPS).any()]
for fn, cols_ in ((f"{P}_groups_selected.aln.fa", showc), (f"{P}_groups_full.aln.fa", fullc)):
    with open(fn, "w") as f:
        for n_, s_ in zip(Cn, Cseq): f.write(f">{n_}\n" + "".join(s_[j] for j in cols_) + "\n")
        for i in order: f.write(f">g{nl['sub'][i]}_i{nl['sub+indel'][i]}|{lab_s[i]}|{cname[i].split('|')[0]}\n" + "".join(cseq[i][j] for j in cols_) + "\n")
out["reading alignments"] = dict(selected_cols_shown=len(showc), full_cols_shown=len(fullc))
json.dump(out, open(f"{P}_v4.json", "w"), indent=1)
print(json.dumps({a: b for a, b in out.items() if not a.startswith("nested") and a != "consensus_pairs_shared_cols_identity"}, indent=1))
print("consensus identity", ci)
for a in ("nested sub", "nested sub+indel"): print(a, {x: y for x, y in out[a].items() if x != "log"})

cols = plt.cm.tab10(np.arange(10)); show = [a for a in ("K2P CpG-masked", "selected sub", "selected sub+indel") if a in pics]
fig, ax = plt.subplots(2, len(show), figsize=(7 * len(show), 12), squeeze=False)
for a_, nm in zip(ax.T, show):
    U, hd = pics[nm]
    for i, s in enumerate(subs): a_[0].scatter(*U[lab == i].T, s=8, color=cols[i % 10], label=s, lw=0, alpha=.8)
    a_[0].set_title(f"UMAP of {nm}\ncoloured by published label; label silhouette on the distance {out[nm]['silh_labels_on_D']:.2f}", fontsize=10)
    cl2 = nl["sub+indel"] if nm == "selected sub+indel" else (nl["sub"] if nm == "selected sub" else hd); c20 = plt.cm.tab20(np.arange(20))
    for q in sorted(set(cl2)): a_[1].scatter(*U[cl2 == q].T, s=8, color="#bbbbbb" if q == -1 else c20[q % 20], lw=0, label=f"{q}")
    a_[1].set_title(("HDBSCAN clusters (number not given)" if nm == "K2P CpG-masked" else "nested groups") + f" on the same picture: {len(set(cl2) - {-1})}", fontsize=10)
    for b in a_: b.set_xticks([]); b.set_yticks([])
    a_[0].legend(markerscale=2, fontsize=7, title="published"); a_[1].legend(markerscale=2, fontsize=7, ncol=2)
fig.suptitle(f"{P}: {N} copies, inputs compared", fontsize=12); plt.tight_layout(); plt.savefig(f"{P}_v4.png", dpi=110)
