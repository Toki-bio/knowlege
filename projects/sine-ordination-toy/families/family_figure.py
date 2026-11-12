"""Presentation figure for one family: copies coloured by their assigned subfamily.

Panels: UMAP of K2P over all columns (CpG masked) | UMAP and t-SNE of a Hamming distance on the columns where the
subfamily consensuses differ (substitutions and indels as characters, pairwise deletion).
Under each panel, measured on the distance and not on the picture: share of copies whose 15 nearest copies are
mostly of their own assigned subfamily, and the silhouette of the assigned subfamilies.
Input : MAFFT alignment of consensuses (rows without '|') and copy bodies (rows 'coord|label').
Output: <P>_figure.png, <P>_figure.json
Usage : python family_figure.py ALN PREFIX "Title"
"""
import sys, json, numpy as np
sys.path.insert(0, ".."); sys.path.insert(0, "/home/toki/ordination_toy")
aln, P, title = sys.argv[1], sys.argv[2], sys.argv[3]; sys.argv = sys.argv[:1]
import sine_inputs_v4 as s4
import ordination_toy_v3 as v3
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

SEED = 1; names, seqs = [], []
for line in open(aln):
    line = line.strip()
    if line.startswith(">"): names.append(line[1:]); seqs.append([])
    elif line: seqs[-1].append(line.upper())
code = {"A": 0, "C": 1, "G": 2, "T": 3}; A = np.array([[code.get(c, -1) for c in "".join(s)] for s in seqs])
iscopy = np.array(["|" in n for n in names]); X = A[iscopy]; Cn = [n for n, c in zip(names, iscopy) if not c]
C = A[~iscopy].copy(); C[C < 0] = s4.GAPS
lab_s = np.array([n.split("|")[1] for n, c in zip(names, iscopy) if c]); subs = [c for c in Cn if c in set(lab_s)]; lab = np.array([subs.index(s) for s in lab_s])
N, L = X.shape; S = s4.code_states(X)
cover = (X >= 0).mean(0); keep = cover >= 0.5
cons = np.array([np.bincount(c[c >= 0], minlength=4).argmax() if (c >= 0).any() else -1 for c in X.T])
cg = np.zeros(L, bool); k = np.where(keep)[0]; ck = cons[k]; cgk = (ck[:-1] == 1) & (ck[1:] == 2); cg[k[:-1][cgk]] = True; cg[k[1:][cgk]] = True
sa, Ua, Ta, _ = s4.evaluate(v3.k2p_pd(X[:, keep & ~cg]), lab, SEED, len(subs))
ft = s4.build_features(S, indel=True); cvar = np.array([len(set(C[:, j].tolist())) > 1 for j in range(L)]); cd = cvar[ft["col"]]
sd, Ud, Td, _ = s4.evaluate(s4.wham(ft["Fm"][:, cd], np.ones(int(cd.sum()))), lab, SEED, len(subs))
out = dict(copies=int(N), subfamilies=subs, all_columns=dict(cols=int((keep & ~cg).sum()), **sa),
           consensus_difference_columns=dict(characters=int(cd.sum()), cols=int(len(set(ft["col"][cd].tolist()))), **sd))
json.dump(out, open(f"{P}_figure.json", "w"), indent=1); print(json.dumps(out, indent=1))
cols = plt.cm.tab10(np.arange(10)); fig, ax = plt.subplots(1, 3, figsize=(19, 6.6))
for a, Z, t, sc in ((ax[0], Ua, f"UMAP, all columns ({out['all_columns']['cols']}), K2P", sa),
                    (ax[1], Ud, f"UMAP, columns where the consensuses differ ({out['consensus_difference_columns']['cols']})", sd),
                    (ax[2], Td, "t-SNE, same columns", sd)):
    for i, s in enumerate(subs): a.scatter(*Z[lab == i].T, s=9, color=cols[i % 10], label=s, lw=0, alpha=.85)
    a.set_title(t, fontsize=11); a.set_xticks([]); a.set_yticks([])
    a.set_xlabel(f"nearest copies mostly of own subfamily: {100 * sc['knn_on_D']:.0f} % of copies;  silhouette {sc['silh_labels_on_D']:.2f}", fontsize=10)
ax[0].legend(markerscale=2, fontsize=8, title="assigned subfamily")
fig.suptitle(f"{title}: {N} copies", fontsize=13); plt.tight_layout(); plt.savefig(f"{P}_figure.png", dpi=110)
