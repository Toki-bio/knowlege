"""Real data, step 4.3: do the published Scalopus aquaticus (saq) Tal subfamilies s1-s9 appear as groups?

Input : saq_all.aln.fa = MAFFT alignment of 900 random copies (100 per subfamily, from Tal repo
        saq/alignments/saq_s*_rand100.aln.fa, ungapped) + 9 subfamily consensuses (rows named s*).
Method (from the toy, page sections 2.6-2.7):
  K2P with pairwise deletion; CpG columns of the majority consensus of the copies masked;
  columns with < 50 % of copies having a base dropped (rare insertions);
  pictures: PCoA, t-SNE (perplexity 40), UMAP (n_neighbors 30, min_dist 0.3) of the distance;
  clusters: Ward k = 9 on each picture (number given) and HDBSCAN (min_cluster_size 20) on UMAP (number not given).
Scores against the published labels: ARI, kNN accuracy (k = 15), silhouette of the published labels on K2P.
Output: saq_ordination.json, saq_ordination.png, saq_hdbscan_vs_labels.tsv, saq_between_copies.tsv
"""
import json, sys, numpy as np
from sklearn.manifold import TSNE
from sklearn.metrics import adjusted_rand_score, silhouette_score
from sklearn.cluster import AgglomerativeClustering, HDBSCAN
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import cross_val_score
import umap, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

SEED = 1; MINOV = 50; CLIP = -0.5 * np.log(1e-6)
names, seqs = [], []
for line in open(sys.argv[1] if len(sys.argv) > 1 else "saq_all.aln.fa"):
    line = line.strip()
    if line.startswith(">"): names.append(line[1:]); seqs.append([])
    elif line: seqs[-1].append(line.upper())
seqs = ["".join(s) for s in seqs]
code = {"A": 0, "C": 1, "G": 2, "T": 3}
A = np.array([[code.get(c, -1) for c in s] for s in seqs])
iscopy = np.array(["|" in n for n in names])
X = A[iscopy]; cname = [n for n, c in zip(names, iscopy) if c]
lab_s = np.array([n.split("|")[1] for n in cname]); subs = sorted(set(lab_s)); lab = np.array([subs.index(s) for s in lab_s])
N, Lc = X.shape
cover = (X >= 0).mean(0); keep = cover >= 0.5
cons = np.array([np.bincount(c[c >= 0], minlength=4).argmax() if (c >= 0).any() else -1 for c in X.T])
cg = np.zeros(Lc, bool); k = np.where(keep)[0]; ck = cons[k]
cgk = (ck[:-1] == 1) & (ck[1:] == 2); cg[k[:-1][cgk]] = True; cg[k[1:][cgk]] = True
use = keep & ~cg
Xu = X[:, use]
info = dict(copies=int(N), subfamilies=subs, align_cols=int(Lc), cols_kept=int(keep.sum()), cpg_cols_masked=int(cg.sum()),
            cols_used=int(use.sum()), gap_frac_used=float((Xu < 0).mean()))

def k2p_pd(X):
    n = len(X); D = np.zeros((n, n)); B = X >= 0
    for i in range(n):
        both = B[i] & B; m = both.sum(1).astype(float)
        diff = (X[i] != X) & both; ts = ((X[i] ^ 2) == X) & both
        with np.errstate(divide="ignore", invalid="ignore"):
            P = ts.sum(1) / m; Q = (diff & ~ts).sum(1) / m; a = 1 - 2 * P - Q; b = 1 - 2 * Q
            D[i] = np.where((a > 1e-6) & (b > 1e-6) & (m >= MINOV),
                            -0.5 * np.log(np.clip(a, 1e-6, None)) - 0.25 * np.log(np.clip(b, 1e-6, None)), CLIP)
    np.fill_diagonal(D, 0); return D
D = k2p_pd(Xu)
sat = D >= CLIP - 1e-9
info.update(saturated_pairs=int(sat.sum() // 2), copies_with_saturated=int(sat.any(1).sum()),
            median_k2p=float(np.median(D[np.triu_indices(N, 1)])),
            silh_published_labels_on_K2P=float(silhouette_score(D, lab, metric="precomputed")))
# per-subfamily mean within vs between K2P
W = {s: float(D[np.ix_(lab == i, lab == i)][np.triu_indices((lab == i).sum(), 1)].mean()) for i, s in enumerate(subs)}
info["mean_k2p_within"] = W

Jc = np.eye(N) - 1 / N; Bm = -0.5 * Jc @ (D ** 2) @ Jc; w, V = np.linalg.eigh(Bm); o = np.argsort(w)[::-1]; w, V = w[o], V[:, o]
info["pcoa_neg_share"] = float(-w[w < 0].sum() / np.abs(w).sum()); info["pcoa12_of_pos"] = float(w[:2].sum() / w[w > 0].sum())
pics = {"PCoA": V[:, :2] * np.sqrt(w[:2]),
        "t-SNE": TSNE(2, metric="precomputed", init="random", perplexity=40, random_state=SEED).fit_transform(D),
        "UMAP": umap.UMAP(n_neighbors=30, min_dist=0.3, metric="precomputed", random_state=SEED).fit_transform(D)}
res = {}
for k_, Z in pics.items():
    cl = AgglomerativeClustering(len(subs), linkage="ward").fit_predict(Z)
    res[k_] = dict(ari_ward=adjusted_rand_score(lab, cl), knn=float(cross_val_score(KNeighborsClassifier(15), Z, lab, cv=5).mean()))
hd = HDBSCAN(min_cluster_size=20).fit_predict(pics["UMAP"])
res["HDBSCAN_on_UMAP"] = dict(ari=adjusted_rand_score(lab, hd), n_clusters=int(len(set(hd) - {-1})), frac_unclustered=float((hd == -1).mean()))
info["results"] = res
# crosstab HDBSCAN cluster x published subfamily
cls = sorted(set(hd))
with open("saq_hdbscan_vs_labels.tsv", "w") as f:
    f.write("hdbscan_cluster\t" + "\t".join(subs) + "\ttotal\n")
    for c in cls:
        row = [int(((hd == c) & (lab == i)).sum()) for i in range(len(subs))]
        f.write(f"{c}\t" + "\t".join(map(str, row)) + f"\t{sum(row)}\n")
# copies whose 15 nearest neighbours in K2P are mostly another subfamily
nn = np.argsort(D, 1)[:, 1:16]
own = (lab[nn] == lab[:, None]).mean(1)
maj = np.array([np.bincount(lab[nn[i]], minlength=len(subs)).argmax() for i in range(N)])
with open("saq_between_copies.tsv", "w") as f:
    f.write("copy\tpublished\tshare_of_15nn_same\tmajority_of_15nn\tmean_k2p_to_15nn\n")
    for i in np.argsort(own):
        if own[i] < 0.5:
            f.write(f"{cname[i]}\t{lab_s[i]}\t{own[i]:.2f}\t{subs[maj[i]]}\t{D[i, nn[i]].mean():.3f}\n")
info["copies_with_15nn_mostly_other_subfamily"] = int((own < 0.5).sum())
info["per_subfamily_share_15nn_same"] = {s: float(own[lab == i].mean()) for i, s in enumerate(subs)}
json.dump(info, open("saq_ordination.json", "w"), indent=1)
print(json.dumps(info, indent=1))

cols = plt.cm.tab10(np.arange(len(subs)))
fig, ax = plt.subplots(1, 4, figsize=(24, 6))
for a, (k_, Z) in zip(ax[:3], pics.items()):
    for i, s in enumerate(subs): a.scatter(*Z[lab == i].T, s=8, color=cols[i], label=s, lw=0, alpha=.8)
    t = f"{k_} of masked K2P\nWard k={len(subs)} ARI {res[k_]['ari_ward']:.2f}, kNN acc {res[k_]['knn']:.2f}"
    if k_ == "PCoA": t += f"\naxes 1+2 = {100*info['pcoa12_of_pos']:.1f}% of +eig; neg share {100*info['pcoa_neg_share']:.0f}%"
    a.set_title(t, fontsize=10); a.set_xticks([]); a.set_yticks([])
ax[0].legend(markerscale=2, fontsize=8, title="published")
Z = pics["UMAP"]; c10 = plt.cm.tab20(np.arange(20))
for j, c in enumerate(cls):
    ax[3].scatter(*Z[hd == c].T, s=8, color="#bbbbbb" if c == -1 else c10[j % 20], lw=0, label="unclustered" if c == -1 else f"c{c}")
ax[3].set_title(f"UMAP, HDBSCAN clusters (number not given)\n{res['HDBSCAN_on_UMAP']['n_clusters']} clusters, ARI vs published {res['HDBSCAN_on_UMAP']['ari']:.2f}", fontsize=10)
ax[3].set_xticks([]); ax[3].set_yticks([]); ax[3].legend(markerscale=2, fontsize=7, ncol=2)
fig.suptitle(f"Scalopus aquaticus Tal SINE: {N} random copies (100 per published subfamily), {info['cols_used']} columns, CpG masked", fontsize=12)
plt.tight_layout(); plt.savefig("saq_ordination.png", dpi=130)
