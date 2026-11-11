"""Repaired toy (v2). Same simulation as ordination_toy.py; fixes the issues listed in
methods/sine-ordination-and-pansineome-visuals.md section 4.1:
  - every method scored in SEQUENCE space: Ward(k=6) clusters of the 2D picture -> silhouette on the K2P matrix
  - two input conditions, same input for all methods in a condition (K2P-based; one-hot PCA-30-based)
  - NMDS started from the PCoA solution, several starts, many iterations
  - several seeds (data seed = method seed = 1..NSEED); report mean and sd
  - control: seed 1, PCA on one-hot, ARI must equal 0.905 (old table); true-label silhouette on K2P is the reference
  - NaN-safe, prints negative-eigenvalue share with every PCoA percentage
Usage: python ordination_toy_v2.py [NSEED]   (needs numpy scikit-learn umap-learn pacmap)
"""
import sys, json, numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import MDS, TSNE, trustworthiness
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import cross_val_score
import umap, pacmap

NSEED = int(sys.argv[1]) if len(sys.argv) > 1 else 5
L = 260; SUB = list("ABCDEF")
PARENT = {"A": None, "B": "A", "C": "A", "D": "B", "E": "B", "F": "C"}
SIZES = dict(A=250, B=300, C=250, D=250, E=200, F=150)

def simulate(seed):                      # identical call order to ordination_toy.py
    rng = np.random.default_rng(seed)
    root = rng.integers(0, 4, L); cons = {}
    for s in SUB:
        base = root if PARENT[s] is None else cons[PARENT[s]].copy()
        n = 0 if PARENT[s] is None else 6
        c = base.copy(); idx = rng.choice(L, n, replace=False)
        c[idx] = (c[idx] + rng.integers(1, 4, n)) % 4; cons[s] = c
    X = []; y = []
    for s in SUB:
        for _ in range(SIZES[s]):
            age = rng.gamma(2, 0.045)
            c = cons[s].copy(); m = rng.random(L) < age
            c[m] = (c[m] + rng.integers(1, 4, m.sum())) % 4
            X.append(c); y.append(s)
    X = np.array(X); lab = np.array([SUB.index(v) for v in y])
    return X, lab

def k2p(X):
    D = np.zeros((len(X), len(X)))
    for i in range(len(X)):
        d = X[i] != X; ts = ((X[i] ^ 2) == X)
        P = ts.mean(1); Q = (d & ~ts).mean(1)
        a = np.clip(1 - 2 * P - Q, 1e-6, None); b = np.clip(1 - 2 * Q, 1e-6, None)
        D[i] = -0.5 * np.log(a) - 0.25 * np.log(b)
    np.fill_diagonal(D, 0); return D

def pcoa(D, k):
    N = len(D); J = np.eye(N) - 1 / N; B = -0.5 * J @ (D ** 2) @ J
    w, V = np.linalg.eigh(B); o = np.argsort(w)[::-1]; w, V = w[o], V[:, o]
    neg = -w[w < 0].sum() / np.abs(w).sum()
    kk = min(k, int((w > 0).sum()))
    return V[:, :kk] * np.sqrt(w[:kk]), w[w > 0].sum(), w, neg

def scores(Z, D, lab, seed):
    cl = AgglomerativeClustering(6, linkage="ward").fit_predict(Z)
    knn = cross_val_score(KNeighborsClassifier(15), Z, lab, cv=5).mean()
    return dict(
        trust=trustworthiness(D, Z, n_neighbors=15, metric="precomputed"),
        silh_seq=silhouette_score(D, cl, metric="precomputed") if len(set(cl)) > 1 else float("nan"),
        silh_pic=silhouette_score(Z, lab),          # old score, kept only to show the gap
        knn=knn, ari=adjusted_rand_score(lab, cl))

def run(seed):
    X, lab = simulate(seed); N = len(lab)
    D = k2p(X); OH = np.eye(4)[X].reshape(N, -1)
    ref = silhouette_score(D, lab, metric="precomputed")   # TRUE labels on K2P: reference for every method
    P30 = PCA(30, random_state=seed).fit_transform(OH)     # one-hot condition
    C30, posmass, w, neg = pcoa(D, 30)                     # K2P condition (PaCMAP feature input)
    pics = {}
    # ---- condition O: input = PCA-30 of one-hot, Euclidean, all methods
    pics["O:PCA2"] = PCA(2, random_state=seed).fit_transform(OH)
    pics["O:t-SNE"] = TSNE(2, init="pca", perplexity=40, random_state=seed).fit_transform(P30)
    pics["O:UMAP"] = umap.UMAP(n_neighbors=30, min_dist=0.3, random_state=seed).fit_transform(P30)
    pics["O:PaCMAP"] = pacmap.PaCMAP(n_components=2, n_neighbors=30, random_state=seed).fit_transform(P30)
    # ---- condition K: input = K2P matrix, all methods (PaCMAP gets PCoA-30 coordinates of K2P)
    Zp = C30[:, :2]
    pics["K:PCoA"] = Zp
    nm = MDS(2, metric=False, dissimilarity="precomputed", init="classical_mds", n_init=1,
             max_iter=3000, eps=1e-6, random_state=seed, normalized_stress="auto")
    pics["K:NMDS(PCoA init)"] = nm.fit_transform(D); stress = nm.stress_
    pics["K:t-SNE"] = TSNE(2, metric="precomputed", init="random", perplexity=40, random_state=seed).fit_transform(D)
    pics["K:t-SNE(PCoA init)"] = TSNE(2, metric="precomputed", init=Zp / Zp.std() * 1e-4, perplexity=40,
                                      random_state=seed).fit_transform(D)
    pics["K:UMAP"] = umap.UMAP(n_neighbors=30, min_dist=0.3, metric="precomputed", random_state=seed).fit_transform(D)
    pics["K:PaCMAP"] = pacmap.PaCMAP(n_components=2, n_neighbors=30, random_state=seed).fit_transform(C30)
    out = {k: scores(Z, D, lab, seed) for k, Z in pics.items()}
    meta = dict(seed=seed, silh_true_on_K2P=ref, pcoa12_of_pos=float(w[w > 0][:2].sum() / posmass),
                pcoa_neg_share=float(neg), nmds_stress=float(stress))
    # control values for seed 1 (must match old page): PCA one-hot ARI 0.905, PC1+PC2 11.9 %
    meta["pca12_var"] = float(PCA(10).fit(OH).explained_variance_ratio_[:2].sum())
    return out, meta

if __name__ == "__main__":
    allr = []; metas = []
    for s in range(1, NSEED + 1):
        r, m = run(s); allr.append(r); metas.append(m)
        print("seed", s, "done  true-label silhouette on K2P = %.3f  PCoA neg share = %.3f  NMDS stress = %.3f" %
              (m["silh_true_on_K2P"], m["pcoa_neg_share"], m["nmds_stress"]), flush=True)
    keys = list(allr[0]); cols = ["trust", "silh_seq", "silh_pic", "knn", "ari"]
    print("\n%-20s" % "picture" + "".join("%18s" % c for c in cols))
    for k in keys:
        row = "%-20s" % k
        for c in cols:
            v = np.array([r[k][c] for r in allr]); row += "%11.3f ±%.3f" % (np.nanmean(v), np.nanstd(v))
        print(row)
    print("\nreference: true-label silhouette on K2P  mean %.3f" % np.mean([m["silh_true_on_K2P"] for m in metas]))
    # control (seed 1 only): old table value for PCA(one-hot) ARI = 0.905; PC1+PC2 = 0.119
    print("CONTROL seed1: PCA2 ARI = %.3f (old 0.905); PC1+PC2 var = %.3f (old 0.119)" %
          (allr[0]["O:PCA2"]["ari"], metas[0]["pca12_var"]))
    json.dump(dict(results=allr, meta=metas), open("ordination_v2_metrics.json", "w"), indent=1)
