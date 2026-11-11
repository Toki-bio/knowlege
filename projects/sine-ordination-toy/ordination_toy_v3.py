"""Harder toy (v3): SINE copies with transitions/CpG decay, indels and 5' truncation.

Same subfamily tree, sizes and per-copy divergence draw as v1/v2; added, with sources:
  - transition/transversion: P(transition) = 2/3 per substitution (my choice, kappa ~2)
  - CpG decay: substitution weight x6 at CpG sites of the copy's consensus, as C->T / G->A
    (Xing et al. 2004 J Mol Biol 344:675: CpG vs non-CpG ratio ~6 in AluY/AluS)
  - indels: 0.15 events per substitution (between Lunter et al. 2006 ~1:10 and Hubley et al. 2022 0.20),
    insertion:deletion = 0.4:0.6, length geometric mean ~1.7, max 20 (Hubley et al. 2022 NAR Genom Bioinf lqac040)
    Deletions become gaps in alignment columns. Insertions are not represented as columns: with pairwise
    deletion they never enter a distance (simplification).
  - 5' truncation: 30 % of copies, truncated part uniform 0-50 % of length (TEvarSim defaults, Miao & Li 2026)
Distance: K2P with pairwise deletion (columns where both copies have a base); pairs with < MINOV shared
columns or saturated K2P get the clip value. One-hot: gap = all-zero.
Arms: ALL copies, and FILTER = drop copies with any saturated K2P cell (automatic, no labels used).
Usage: python ordination_toy_v3.py [NSEED]    env: OUT=file.json
"""
import os, sys, json, numpy as np
from sklearn.decomposition import PCA
from sklearn.manifold import MDS, TSNE, trustworthiness
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.cluster import AgglomerativeClustering
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import cross_val_score
import umap, pacmap

NSEED = int(sys.argv[1]) if len(sys.argv) > 1 else 5
L = 260; SUB = list("ABCDEF"); GAP = -1
PARENT = {"A": None, "B": "A", "C": "A", "D": "B", "E": "B", "F": "C"}
SIZES = dict(A=250, B=300, C=250, D=250, E=200, F=150)
P_TS = 2 / 3; CPG_W = 6.0; INDEL_PER_SUB = 0.15; P_INS = 0.4; TRUNC_MAX = 0.5
TRUNC_FRAC = float(os.environ.get("TRUNC_FRAC", "0.3"))   # 0 = no 5' truncation (test arm)
MASK_CPG = os.environ.get("MASK_CPG") == "1"              # mask CpG columns of the family majority consensus (COSEG practice)
from sklearn.cluster import HDBSCAN
MINOV = 50; CLIP = -0.5 * np.log(1e-6)   # 0=A 1=C 2=G 3=T ; transition partner = x ^ 2

def indel_len(rng):
    return int(min(20, rng.geometric(0.6)))       # mean ~1.67

def make_copy(cons, d, rng):
    c = cons.copy()
    cpg = np.zeros(L, bool)
    cg = (cons[:-1] == 1) & (cons[1:] == 2); cpg[:-1] |= cg; cpg[1:] |= cg
    w = np.where(cpg, CPG_W, 1.0); w = w / w.mean()
    m = rng.random(L) < np.clip(d * w, 0, 0.95)
    nsub = int(m.sum())
    for i in np.where(m)[0]:
        if cpg[i]:
            c[i] = 3 if cons[i] == 1 else 0          # C->T, G->A
        elif rng.random() < P_TS:
            c[i] = c[i] ^ 2
        else:
            c[i] = (c[i] + rng.choice([1, 3])) % 4   # transversion
    nindel = rng.poisson(INDEL_PER_SUB * nsub)
    for _ in range(nindel):
        if rng.random() >= P_INS:                    # deletion -> gap columns
            p = rng.integers(0, L); c[p:p + indel_len(rng)] = GAP
    if rng.random() < TRUNC_FRAC:                    # 5' truncation
        c[:int(rng.uniform(0, TRUNC_MAX) * L)] = GAP
    return c

def simulate(seed):
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
            X.append(make_copy(cons[s], rng.gamma(2, 0.045), rng)); y.append(SUB.index(s))
    return np.array(X), np.array(y)

def k2p_pd(X):
    N = len(X); D = np.zeros((N, N)); B = X >= 0
    for i in range(N):
        both = B[i] & B; n = both.sum(1).astype(float)
        diff = (X[i] != X) & both; ts = ((X[i] ^ 2) == X) & both
        with np.errstate(divide="ignore", invalid="ignore"):
            P = ts.sum(1) / n; Q = (diff & ~ts).sum(1) / n
            a = 1 - 2 * P - Q; b = 1 - 2 * Q
            d = np.where((a > 1e-6) & (b > 1e-6) & (n >= MINOV), -0.5 * np.log(np.clip(a, 1e-6, None)) - 0.25 * np.log(np.clip(b, 1e-6, None)), CLIP)
        D[i] = d
    np.fill_diagonal(D, 0); return D

def pcoa(D, k):
    N = len(D); J = np.eye(N) - 1 / N; B = -0.5 * J @ (D ** 2) @ J
    w, V = np.linalg.eigh(B); o = np.argsort(w)[::-1]; w, V = w[o], V[:, o]
    kk = min(k, int((w > 0).sum()))
    return V[:, :kk] * np.sqrt(w[:kk]), -w[w < 0].sum() / np.abs(w).sum()

def scores(Z, D, lab):
    cl = AgglomerativeClustering(6, linkage="ward").fit_predict(Z)
    return dict(trust=trustworthiness(D, Z, n_neighbors=15, metric="precomputed"),
                silh_seq=silhouette_score(D, cl, metric="precomputed"),
                knn=cross_val_score(KNeighborsClassifier(15), Z, lab, cv=5).mean(),
                ari=adjusted_rand_score(lab, cl))

def pictures(X, D, seed):
    N = len(X); OH = np.zeros((N, L * 4)); idx = np.where(X >= 0)
    OH[idx[0], idx[1] * 4 + X[idx]] = 1
    P30 = PCA(30, random_state=seed).fit_transform(OH); C30, neg = pcoa(D, 30); Zp = C30[:, :2]
    nm = MDS(2, metric_mds=False, metric="precomputed", init="classical_mds", n_init=1, max_iter=3000, eps=1e-6, random_state=seed)
    pics = {"O:PCA2": PCA(2, random_state=seed).fit_transform(OH),
            "O:UMAP": umap.UMAP(n_neighbors=30, min_dist=0.3, random_state=seed).fit_transform(P30),
            "O:PaCMAP": pacmap.PaCMAP(n_components=2, n_neighbors=30, random_state=seed).fit_transform(P30),
            "K:PCoA": Zp, "K:NMDS": nm.fit_transform(D),
            "K:t-SNE": TSNE(2, metric="precomputed", init="random", perplexity=40, random_state=seed).fit_transform(D),
            "K:UMAP": umap.UMAP(n_neighbors=30, min_dist=0.3, metric="precomputed", random_state=seed).fit_transform(D),
            "K:PaCMAP": pacmap.PaCMAP(n_components=2, n_neighbors=30, random_state=seed).fit_transform(C30)}
    return pics, dict(pcoa_neg_share=float(neg), nmds_stress=float(nm.stress_))

def hdb(Z, lab, metric="euclidean", mcs=20):
    cl = HDBSCAN(min_cluster_size=mcs, metric=metric).fit_predict(Z)
    return dict(ari=adjusted_rand_score(lab, cl), n_clusters=int(len(set(cl) - {-1})), frac_noise=float((cl == -1).mean()))

def run(seed):
    X, lab = simulate(seed)
    if MASK_CPG:
        cons = np.array([np.bincount(c[c >= 0], minlength=4).argmax() for c in X.T])
        cg = (cons[:-1] == 1) & (cons[1:] == 2); m = np.zeros(L, bool); m[:-1] |= cg; m[1:] |= cg
        X = X.copy(); X[:, m] = GAP
    D = k2p_pd(X)
    sat = (D >= CLIP - 1e-9); bad = sat.any(1)
    meta = dict(seed=seed, frac_gap=float((X < 0).mean()), frac_truncated=float((X[:, 0] < 0).mean()),
                n_sat_cells=int(sat.sum() // 2), n_copies_with_sat=int(bad.sum()))
    out = {}
    for arm, keep in (("ALL", np.ones(len(X), bool)), ("FILTER", ~bad)):
        Xa, la, Da = X[keep], lab[keep], D[np.ix_(keep, keep)]
        pics, m = pictures(Xa, Da, seed)
        out[arm] = {k: scores(Z, Da, la) for k, Z in pics.items()}
        if arm == "ALL":   # number of subfamilies NOT given
            OHa = np.zeros((len(Xa), L * 4)); ii = np.where(Xa >= 0); OHa[ii[0], ii[1] * 4 + Xa[ii]] = 1
            meta["HDBSCAN"] = {"K:UMAP picture": hdb(pics["K:UMAP"], la), "K:t-SNE picture": hdb(pics["K:t-SNE"], la),
                               "K2P matrix": hdb(Da.copy(), la, metric="precomputed"),   # HDBSCAN overwrites a precomputed matrix in place
                               "one-hot PCA-100 (Zhou 2019)": hdb(PCA(100, random_state=seed).fit_transform(OHa), la)}
        meta[arm] = dict(n=int(keep.sum()), silh_true_on_K2P=float(silhouette_score(Da, la, metric="precomputed")), **m)
    return out, meta

if __name__ == "__main__":
    allr, metas = [], []
    for s in range(1, NSEED + 1):
        r, m = run(s); allr.append(r); metas.append(m)
        print("seed", s, json.dumps(m), flush=True)
        json.dump(dict(results=allr, meta=metas), open(os.environ.get("OUT", "ordination_v3_metrics.json"), "w"), indent=1)  # checkpoint
    for arm in ("ALL", "FILTER"):
        print(f"\n== {arm} ==   true-label silhouette on K2P mean {np.mean([m[arm]['silh_true_on_K2P'] for m in metas]):.3f}")
        print("%-12s" % "picture" + "".join("%18s" % c for c in ["trust", "silh_seq", "knn", "ari"]))
        for k in allr[0][arm]:
            print("%-12s" % k + "".join("%11.3f ±%.3f" % (np.mean([r[arm][k][c] for r in allr]), np.std([r[arm][k][c] for r in allr]))
                                         for c in ["trust", "silh_seq", "knn", "ari"]))
        print("ARI per seed:", {k: [round(r[arm][k]["ari"], 3) for r in allr] for k in allr[0][arm]})
    print("\n== HDBSCAN (min_cluster_size 20, k not given; noise counted as its own label) ==")
    for k in metas[0]["HDBSCAN"]:
        v = [m["HDBSCAN"][k] for m in metas]
        print("%-30s ARI %.3f ±%.3f  clusters %s  noise %.2f" % (k, np.mean([x["ari"] for x in v]), np.std([x["ari"] for x in v]),
                                                              [x["n_clusters"] for x in v], np.mean([x["frac_noise"] for x in v])))
    json.dump(dict(results=allr, meta=metas), open(os.environ.get("OUT", "ordination_v3_metrics.json"), "w"), indent=1)
