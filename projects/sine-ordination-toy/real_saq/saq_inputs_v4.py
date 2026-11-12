"""saq, 866 copies: inputs other than unweighted K2P (functions and method: ../sine_inputs_v4.py).

Arms: K2P CpG-masked (baseline of 2.8); all minor-allele features; features selected by covariation between
positions (label-free), substitutions only and with indels as characters; the same weighted; PCA with parallel
analysis; nested re-selection inside each cluster. One flagged arm uses the columns where the 9 published
consensuses differ: it is NOT independent of the labels (copies were assigned to those consensuses).
Checks: (1) do the label-free selected positions coincide with positions where the published consensuses differ;
(2) split-half: are the same positions selected in two disjoint random halves of the copies.
Output: saq_inputs_v4.json, saq_inputs_v4_features.tsv, saq_inputs_v4_crosstabs.txt, saq_inputs_v4.png
"""
import sys, json, numpy as np
sys.path.insert(0, "..")
import sine_inputs_v4 as s4
import ordination_toy_v3 as v3
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

SEED = 1
names, seqs = [], []
for line in open("saq_all.aln.fa"):
    line = line.strip()
    if line.startswith(">"): names.append(line[1:]); seqs.append([])
    elif line: seqs[-1].append(line.upper())
code = {"A": 0, "C": 1, "G": 2, "T": 3}; A = np.array([[code.get(c, -1) for c in "".join(s)] for s in seqs])
iscopy = np.array(["|" in n for n in names]); X = A[iscopy]; cname = [n for n, c in zip(names, iscopy) if c]
Cn = [n for n, c in zip(names, iscopy) if not c]; C = A[~iscopy].copy(); C[C < 0] = s4.GAPS
lab_s = np.array([n.split("|")[1] for n in cname]); subs = sorted(set(lab_s)); lab = np.array([subs.index(s) for s in lab_s])
N, L = X.shape; S = s4.code_states(X); LET = "ACGT-"
out = dict(copies=int(N), align_cols=int(L), subfamilies=subs); xt = open("saq_inputs_v4_crosstabs.txt", "w")

def crosstab(title, cl):
    xt.write(f"## {title}\ncluster\t" + "\t".join(subs) + "\ttotal\n")
    for c in sorted(set(cl)):
        row = [int(((cl == c) & (lab == i)).sum()) for i in range(len(subs))]; xt.write(f"{c}\t" + "\t".join(map(str, row)) + f"\t{sum(row)}\n")
    xt.write("\n")

pics = {}
def arm(name, D):
    sc, U, T, hd = s4.evaluate(D, lab, SEED, len(subs)); out[name] = sc; pics[name] = (U, hd); crosstab(f"{name}: HDBSCAN on UMAP x published", hd)
    print(name, json.dumps(sc), flush=True)

# baseline exactly as in saq_ordination.py
cover = (X >= 0).mean(0); keep = cover >= 0.5
cons = np.array([np.bincount(c[c >= 0], minlength=4).argmax() if (c >= 0).any() else -1 for c in X.T])
cg = np.zeros(L, bool); k = np.where(keep)[0]; ck = cons[k]; cgk = (ck[:-1] == 1) & (ck[1:] == 2); cg[k[:-1][cgk]] = True; cg[k[1:][cgk]] = True
arm("K2P CpG-masked", v3.k2p_pd(X[:, keep & ~cg]))

res = {}
for tag, indel in (("sub", False), ("sub+indel", True)):
    ft = s4.build_features(S, indel=indel); w, c = s4.cov_weights(ft, S, SEED); sel = w > 0; res[tag] = (ft, w, sel)
    out[f"selection {tag}"] = dict(features=int(ft["Fm"].shape[1]), selected=int(sel.sum()), selected_indel=int((sel & (ft["kind"] != "sub")).sum()),
                                   threshold_z2=c, selected_on_CpG_cols=int(cg[ft["col"][sel]].sum()))
    if tag == "sub": arm("feat-all", s4.wham(ft["Fm"], np.ones(ft["Fm"].shape[1])))
    if sel.sum() >= 2:
        arm(f"cov-selected {tag}", s4.wham(ft["Fm"][:, sel], np.ones(int(sel.sum()))))
        arm(f"cov-weighted {tag}", s4.wham(ft["Fm"][:, sel], np.sqrt(w[sel])))
        for sh in (0.5, 0.75, 0.9): arm(f"heavy{int(sh * 100)} {tag}", s4.wham(ft["Fm"], s4.heavy(sel, sh)))
ft, w, sel = res["sub+indel"]
Z, kk = s4.pca_pa(res["sub"][0], SEED); out["PCA-pa axes"] = kk
if kk: arm("PCA-pa", s4.euclid(Z))

# consensus comparison: is the selected feature a state in which the published consensuses differ?
def cons_info(j, kind, state):
    cs = C[:, j]; var = len(set(cs.tolist())) > 1
    if kind == "ins": hit = var and (cs != s4.GAPS).any() and (cs == s4.GAPS).any()
    else: hit = var and (cs == state).any()
    return var, bool(hit), "".join(LET[x] for x in cs)
cvar = np.array([len(set(C[:, j].tolist())) > 1 for j in range(L)])
fcols = set(int(j) for j in ft["col"])
with open("saq_inputs_v4_features.tsv", "w") as f:
    f.write("column\tkind\tstate\tfreq\tweight\tselected\tCpG_col\tconsensuses_differ_here\tstate_is_a_consensus_state\tconsensus_states_" + ",".join(Cn) + "\t" + "\t".join("freq_" + s for s in subs) + "\n")
    nsel_hit = 0
    for i in np.argsort(-w):
        j = int(ft["col"][i]); var, hit, cs = cons_info(j, ft["kind"][i], ft["state"][i]); nsel_hit += bool(sel[i] and hit)
        fr = [np.nanmean(ft["Fm"][lab == q, i]) for q in range(len(subs))]
        f.write(f"{j}\t{ft['kind'][i]}\t{LET[ft['state'][i]] if ft['state'][i] >= 0 else 'base'}\t{np.nanmean(ft['Fm'][:, i]):.3f}\t{w[i]:.1f}\t{int(sel[i])}\t{int(cg[j])}\t{int(var)}\t{int(hit)}\t{cs}\t" + "\t".join(f"{x:.2f}" for x in fr) + "\n")
selcols = set(int(j) for j in ft["col"][sel]); cv_feat = set(j for j in fcols if cvar[j])
out["consensus check (sub+indel)"] = dict(selected_features=int(sel.sum()), selected_features_that_are_a_consensus_difference=int(nsel_hit),
                                         columns_where_consensuses_differ=int(cvar.sum()), of_those_with_a_feature=len(cv_feat),
                                         of_those_selected=len(cv_feat & selcols), selected_columns=len(selcols), selected_columns_where_consensuses_do_not_differ=len(selcols - set(np.where(cvar)[0].tolist())))

# flagged arm: columns where the published consensuses differ (not independent of the labels)
cd = np.array([cvar[j] for j in ft["col"]])
if cd.sum() >= 2: arm("consensus-difference columns (NOT label-free)", s4.wham(ft["Fm"][:, cd], np.ones(int(cd.sum()))))

# split-half replication of the selection
rng = np.random.default_rng(7); halves = []
for rep in range(5):
    p = rng.permutation(N); h = []
    for part in (p[:N // 2], p[N // 2:]):
        fh = s4.build_features(S[part], indel=True); wh, _ = s4.cov_weights(fh, S[part], SEED + rep)
        h.append(set((int(a), str(b), int(c_)) for a, b, c_ in zip(fh["col"][wh > 0], fh["kind"][wh > 0], fh["state"][wh > 0])))
    halves.append(dict(a=len(h[0]), b=len(h[1]), both=len(h[0] & h[1])))
out["split-half selection (5 random splits)"] = halves

# nested
for tag, indel in (("sub", False), ("sub+indel", True)):
    log = []; lv = s4.nested(S, SEED, indel=indel, log=log); nl = s4.leaves_to_labels(lv, N)
    from sklearn.metrics import adjusted_rand_score
    out[f"nested {tag}"] = dict(leaves=len(lv), sizes=[int(len(v)) for _, v in sorted(lv.items())], ari_vs_published=float(adjusted_rand_score(lab, nl)), log=log)
    crosstab(f"nested {tag}: leaves x published (leaf order = {sorted(lv)})", nl)
    if indel: nested_lab = nl
xt.close(); json.dump(out, open("saq_inputs_v4.json", "w"), indent=1)
print(json.dumps({a: b for a, b in out.items() if not a.startswith("nested")}, indent=1)); print({a: {x: y for x, y in b.items() if x != "log"} for a, b in out.items() if a.startswith("nested")})

cols = plt.cm.tab10(np.arange(len(subs))); show = [a for a in ("K2P CpG-masked", "heavy75 sub", "heavy75 sub+indel", "cov-selected sub+indel") if a in pics]
fig, ax = plt.subplots(2, len(show), figsize=(7 * len(show), 12), squeeze=False)
for a_, nm in zip(ax.T, show):
    U, hd = pics[nm]
    for i, s in enumerate(subs): a_[0].scatter(*U[lab == i].T, s=8, color=cols[i], label=s, lw=0, alpha=.8)
    a_[0].set_title(f"UMAP of {nm}\ncoloured by published label; label silhouette on the distance {out[nm]['silh_labels_on_D']:.2f}", fontsize=10)
    cl2 = nested_lab if nm == "cov-selected sub+indel" else hd; c20 = plt.cm.tab20(np.arange(20))
    for q in sorted(set(cl2)): a_[1].scatter(*U[cl2 == q].T, s=8, color="#bbbbbb" if q == -1 else c20[q % 20], lw=0, label=f"{q}")
    a_[1].set_title(("nested leaves" if nm == "cov-selected sub+indel" else "HDBSCAN clusters (number not given)") + f" on the same picture: {len(set(cl2) - {-1})}", fontsize=10)
    for b in a_: b.set_xticks([]); b.set_yticks([])
    a_[0].legend(markerscale=2, fontsize=7, title="published"); a_[1].legend(markerscale=2, fontsize=7, ncol=2)
fig.suptitle(f"Scalopus aquaticus Tal SINE, {N} copies: input compared", fontsize=12); plt.tight_layout(); plt.savefig("saq_inputs_v4.png", dpi=120)
