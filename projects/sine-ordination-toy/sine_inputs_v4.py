"""Inputs other than unweighted K2P for pictures of SINE copies (v4).

Question: can a different input show subfamilies that an unweighted distance over all columns cannot?
Idea: a subfamily is defined by a few diagnostic positions. Those positions CO-VARY across copies (a copy that has
one usually has the others); private mutations do not. So positions can be selected / weighted without labels.

Inputs compared (same downstream for all):
  K2P            K2P, pairwise deletion (v3), with and without CpG mask
  feat-all       all minor-allele features (frequency >= MIN_FREQ), unweighted Hamming, pairwise deletion
  cov-selected   only features with significant covariation with a feature at another position; unweighted
  cov-weighted   same, weight = sqrt(summed excess z^2)
  heavy50/75/90  ALL features; the selected ones together carry 50 / 75 / 90 % of the weight
  cov-sel+indel  cov-selected, with internal gaps / insertions as characters
  PCA-pa         PCA of the mean-imputed centred feature matrix, number of axes by parallel analysis
  oracle         (toy only) Hamming on the true diagnostic columns = upper bound
  nested         cov-selected -> UMAP -> HDBSCAN proposes parts -> parts merged until every two differ at >= 2
                 positions (frequency difference >= 0.5) -> features re-selected INSIDE each part, recursively

Covariation: phi correlation between two features over copies where both are observed, after subtracting from each
feature its mean among copies of similar load (older copies carry more of everything, which alone makes all
noise sites co-vary). Threshold = largest z^2 in data permuted within load bins (3 permutations).
Pairs closer than 3 columns (25 if an indel is involved) are not counted.

Toy usage: python sine_inputs_v4.py NSEED    env: DIAG=6|3|2  NOSUB=1  TRUNC_FRAC  OUT=file.json
"""
import os, sys, json, numpy as np
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.cluster import AgglomerativeClustering, HDBSCAN
import umap

MISS, GAPS = -2, 4
MIN_FREQ = float(os.environ.get("MIN_FREQ", "0.03")); MIN_N = 10; MIN_COV = 0.5

def code_states(X):
    """-1 = gap in. Terminal gap runs -> MISS (truncation / not observed); internal gaps -> state 4."""
    S = X.copy(); S[S < 0] = GAPS
    for r in S:
        nz = np.where(r != GAPS)[0]
        if len(nz) == 0: r[:] = MISS; continue
        r[:nz[0]] = MISS; r[nz[-1] + 1:] = MISS
    return S

def build_features(S, indel=False, min_freq=None):
    """Binary features (column, minor state); NaN where the copy is not observed. Also the private-mutation load."""
    min_freq = MIN_FREQ if min_freq is None else min_freq
    N, L = S.shape; vec, col, kind, state = [], [], [], []
    priv = np.zeros(N); seen = np.zeros(N)
    for j in range(L):
        c = S[:, j]; ok = c != MISS; n = int(ok.sum())
        if n < MIN_COV * N: continue
        maj = int(np.bincount(c[ok], minlength=5).argmax())
        cand = [("ins", None)] if maj == GAPS else [("sub", a) for a in range(4) if a != maj] + [("del", GAPS)]
        used = np.zeros(N, bool)
        for kd, a in cand:
            x = ((c != GAPS) if kd == "ins" else (c == a)) & ok
            if x.sum() >= MIN_N and x.sum() / n >= min_freq and (indel or kd == "sub"):
                vec.append(np.where(ok, x.astype(float), np.nan)); col.append(j); kind.append(kd); state.append(-1 if a is None else a)
                used |= x
        if maj != GAPS:
            priv += ok & (c != maj) & (c != GAPS) & ~used; seen += ok
    if not vec: return None
    Fm = np.array(vec).T; col = np.array(col); kind = np.array(kind); state = np.array(state)
    keep = np.ones(len(col), bool); last = {}                       # merge runs of one indel over adjacent columns
    for f in np.argsort(col, kind="stable"):
        if kind[f] == "sub": continue
        p = last.get(kind[f])
        if p is not None and col[f] - p[1] <= 2:
            m = ~np.isnan(Fm[:, f]) & ~np.isnan(Fm[:, p[0]])
            if m.sum() > 20 and np.std(Fm[m, f]) > 0 and np.std(Fm[m, p[0]]) > 0 and np.corrcoef(Fm[m, f], Fm[m, p[0]])[0, 1] > 0.8:
                keep[f] = False; last[kind[f]] = (p[0], col[f]); continue
        last[kind[f]] = (f, col[f])
    return dict(Fm=Fm[:, keep], col=col[keep], kind=kind[keep], state=state[keep], load=priv / np.maximum(seen, 1))

def total_load(S, excl=()):
    """Share of observed columns where the copy differs from the majority base, columns in excl left out."""
    N, L = S.shape; d = np.zeros(N); seen = np.zeros(N); excl = set(excl)
    for j in range(L):
        if j in excl: continue
        c = S[:, j]; ok = c != MISS
        if ok.sum() < MIN_COV * N: continue
        maj = int(np.bincount(c[ok], minlength=5).argmax())
        if maj == GAPS: continue
        d += ok & (c != maj) & (c != GAPS); seen += ok
    return d / np.maximum(seen, 1)

def _bins(load):
    """Bins of copies of similar load: quantile bins, each cut further where it spans more than 0.2 x median load
    (the oldest copies otherwise share one wide bin, and their common CpG decay looks like covariation)."""
    nb = int(min(20, max(2, len(load) // 30))); r = np.argsort(np.argsort(load, kind="stable"), kind="stable")
    q = (r * nb // len(load)).astype(int); wd = 0.2 * max(float(np.median(load)), 1e-9)
    key = q * 100000 + np.floor(load / wd).astype(int); _, b = np.unique(key, return_inverse=True)
    return b, int(b.max()) + 1

def _z2(Fm, load, col, isindel):
    M = ~np.isnan(Fm); R = np.zeros(Fm.shape); b, nb = _bins(load)
    for q in range(nb):                                   # residual = value minus the mean of copies of similar load
        rows = b == q; sub = Fm[rows]; mu = np.nanmean(np.where(np.isnan(sub), np.nan, sub), 0) if rows.any() else 0
        R[rows] = np.nan_to_num(sub - mu)
    Mf = M.astype(float); num = R.T @ R; v1 = (R ** 2).T @ Mf; n = Mf.T @ Mf
    z2 = num ** 2 / (v1 * v1.T + 1e-12) * np.maximum(n - 2, 0)
    dc = np.abs(col[:, None] - col[None, :]); ind = isindel[:, None] | isindel[None, :]
    z2[~(np.where(ind, dc >= 25, dc >= 3) & (n >= 50))] = 0
    return z2

def _permute(Fm, load, rng):
    P = Fm.copy(); b, nb = _bins(load); M = ~np.isnan(Fm)
    for f in range(Fm.shape[1]):
        for q in range(nb):
            idx = np.where(M[:, f] & (b == q))[0]
            if len(idx) > 1: P[idx, f] = Fm[rng.permutation(idx), f]
    return P

def cov_weights(ft, S, seed=0, nperm=3):
    """Summed excess z^2 per feature (0 = no significant partner at another position), and the threshold.
    Two passes: (1) load = private mutations; (2) load = all differences from the majority outside the columns
    selected in pass 1 and their neighbours (a less noisy estimate of copy age, free of the subfamily signal)."""
    Fm, col = ft["Fm"], ft["col"]; isindel = ft["kind"] != "sub"
    if Fm.shape[1] < 2: return np.zeros(Fm.shape[1]), float("nan")
    rng = np.random.default_rng(1000 + seed); load = ft["load"]
    for _ in range(2):
        c = max(_z2(_permute(Fm, load, rng), load, col, isindel).max() for _ in range(nperm))
        z2 = _z2(Fm, load, col, isindel); w = np.maximum(z2 - c, 0).sum(1); ft["sig"] = z2 > c
        load = total_load(S, [int(j) + k for j in col[w > 0] for k in (-1, 0, 1)])
    return w, float(c)

def heavy(sel, share):
    """Weights over ALL features: 1 for the others, and the selected ones together carry `share` of the total weight."""
    w = np.ones(len(sel)); ns = int(sel.sum())
    if 0 < ns < len(sel): w[sel] = share / (1 - share) * (len(sel) - ns) / ns
    return w

def wham(Fm, w):
    """Weighted Hamming distance with pairwise deletion."""
    M = (~np.isnan(Fm)).astype(float); X0 = np.nan_to_num(Fm); XW = X0 * w
    num = XW @ M.T + M @ XW.T - 2 * XW @ X0.T; den = (M * w) @ M.T
    D = np.where(den > 0, num / np.maximum(den, 1e-12), 1.0); D = np.clip((D + D.T) / 2, 0, None); np.fill_diagonal(D, 0)
    return D

def pca_pa(ft, seed=0, nperm=2):
    """PCA of mean-imputed centred features; axes kept = singular values above the largest of permuted data."""
    def sv(F):
        X = F - np.nanmean(F, 0); X = np.nan_to_num(X); return np.linalg.svd(X, full_matrices=False)
    rng = np.random.default_rng(2000 + seed); U, s, _ = sv(ft["Fm"])
    thr = max(sv(_permute(ft["Fm"], ft["load"], rng))[1][0] for _ in range(nperm)); k = int((s > thr).sum())
    return (U[:, :k] * s[:k] if k else None), k

def euclid(Z):
    g = (Z ** 2).sum(1); D = np.sqrt(np.clip(g[:, None] + g[None, :] - 2 * Z @ Z.T, 0, None)); np.fill_diagonal(D, 0); return D

def embed(D, seed):
    n = len(D)
    U = umap.UMAP(n_neighbors=min(30, n - 1), min_dist=0.3, metric="precomputed", random_state=seed).fit_transform(D)
    T = TSNE(2, metric="precomputed", init="random", perplexity=min(40, (n - 1) / 3), random_state=seed).fit_transform(D)
    return U, T

def evaluate(D, lab, seed, k=None):
    """Scores of one input. kNN and silhouette are on the distance itself, not on a picture."""
    n = len(D); k = len(set(lab)) if k is None else k; U, T = embed(D, seed)
    Dn = D + np.diag(np.full(n, np.inf)); nn = np.argsort(Dn + 1e-9 * np.random.default_rng(seed).random(D.shape), 1)[:, :15]
    pred = np.array([np.bincount(lab[r], minlength=lab.max() + 1).argmax() for r in nn])
    hd = HDBSCAN(min_cluster_size=20).fit_predict(U)
    return dict(knn_on_D=float((pred == lab).mean()), silh_labels_on_D=float(silhouette_score(D, lab, metric="precomputed")) if k > 1 else None,
                ari_umap_ward=float(adjusted_rand_score(lab, AgglomerativeClustering(k, linkage="ward").fit_predict(U))) if k > 1 else None,
                ari_tsne_ward=float(adjusted_rand_score(lab, AgglomerativeClustering(k, linkage="ward").fit_predict(T))) if k > 1 else None,
                hdb_ari=float(adjusted_rand_score(lab, hd)), hdb_n=int(len(set(hd) - {-1})), hdb_noise=float((hd == -1).mean())), U, T, hd

DIAG_DF = 0.5; MIN_DIAG = 2

def diag_positions(ft, sel, ia, ib):
    """Selected features whose frequency differs by >= DIAG_DF between two sets of copies AND that co-vary
    significantly, in the whole node, with another such feature (a part made of copies that merely share two
    unrelated private mutations fails this). Positions closer than 3 columns (25 if an indel) count once."""
    with np.errstate(all="ignore"):
        d = np.abs(np.nanmean(ft["Fm"][ia][:, sel], 0) - np.nanmean(ft["Fm"][ib][:, sel], 0)) >= DIAG_DF
    fi = np.where(sel)[0][d]; fi = fi[[bool(ft["sig"][f, fi].any()) for f in fi]] if len(fi) else fi   # must co-vary with another of them in the node
    col = ft["col"][fi]; ind = ft["kind"][fi] != "sub"; o = np.argsort(col); out = []; last = None
    for j, i in zip(col[o], ind[o]):
        if last is None or j - last[0] >= (25 if (i or last[1]) else 3): out.append(int(j))
        last = (j, i)
    return out

def nested(S, seed, indel=False, min_size=60, log=None, idx=None, path="r"):
    """Recursive: select co-varying features inside the set; propose parts by HDBSCAN on UMAP; merge parts until every
    two of them differ at >= MIN_DIAG positions (frequency difference >= DIAG_DF); repeat inside each part."""
    idx = np.arange(len(S)) if idx is None else idx; n = len(idx)
    def leaf(why):
        if log is not None: log.append(dict(node=path, n=int(n), leaf=why))
        return {path: idx}
    if n < min_size: return leaf("small")
    ft = build_features(S[idx], indel=indel)
    if ft is None: return leaf("no features")
    w, c = cov_weights(ft, S[idx], seed); sel = w > 0
    if sel.sum() < 2: return leaf("no co-varying features")
    D = wham(ft["Fm"][:, sel], np.ones(int(sel.sum())))
    U = umap.UMAP(n_neighbors=min(30, n - 1), min_dist=0.3, metric="precomputed", random_state=seed).fit_transform(D)
    hd = HDBSCAN(min_cluster_size=max(20, int(0.03 * n))).fit_predict(U); cl = sorted(set(hd) - {-1})
    if len(cl) < 2: return leaf(f"{int(sel.sum())} co-varying features, one cluster")
    for i in np.where(hd == -1)[0]: hd[i] = cl[int(np.argmin([D[i, hd == q].mean() for q in cl]))]
    parts = [np.where(hd == q)[0] for q in cl]; proposed = len(parts)
    while len(parts) > 1:
        best = min(((len(diag_positions(ft, sel, parts[a], parts[b])), a, b) for a in range(len(parts)) for b in range(a + 1, len(parts))))
        if best[0] >= MIN_DIAG: break
        parts[best[1]] = np.concatenate([parts[best[1]], parts[best[2]]]); parts.pop(best[2])
    if len(parts) < 2: return leaf(f"{int(sel.sum())} co-varying features, {proposed} proposed parts not separated by {MIN_DIAG} positions")
    if log is not None:
        log.append(dict(node=path, n=int(n), features=int(ft["Fm"].shape[1]), selected=int(sel.sum()), proposed=proposed,
                        split=[int(len(p_)) for p_ in parts],
                        diag_positions_between_parts={f"{a}-{b}": diag_positions(ft, sel, parts[a], parts[b]) for a in range(len(parts)) for b in range(a + 1, len(parts))}))
    out = {}
    for q, p_ in enumerate(parts): out.update(nested(S, seed, indel, min_size, log, idx[p_], f"{path}.{q}"))
    return out

def leaves_to_labels(leaves, n):
    lab = np.zeros(n, int)
    for i, (k, v) in enumerate(sorted(leaves.items())): lab[v] = i
    return lab

# ----------------------------------------------------------------------------------------------- toy
if __name__ == "__main__":
    import ordination_toy_v3 as v3
    NSEED = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    DIAG = int(os.environ.get("DIAG", "6")); NOSUB = os.environ.get("NOSUB") == "1"; OUT = os.environ.get("OUT", "inputs_v4.json")

    def simulate(seed):
        rng = np.random.default_rng(seed); root = rng.integers(0, 4, v3.L); cons = {}; dcols = set()
        for s in v3.SUB:
            base = root if v3.PARENT[s] is None else cons[v3.PARENT[s]].copy()
            n = 0 if (v3.PARENT[s] is None or NOSUB) else DIAG
            c = base.copy(); idx = rng.choice(v3.L, n, replace=False); c[idx] = (c[idx] + rng.integers(1, 4, n)) % 4
            cons[s] = c; dcols |= set(int(i) for i in idx)
        X, y = [], []
        for s in v3.SUB:
            for _ in range(v3.SIZES[s]): X.append(v3.make_copy(cons[s], rng.gamma(2, 0.045), rng)); y.append(v3.SUB.index(s))
        return np.array(X), np.array(y), sorted(dcols)

    def cpg_mask(X):
        cons = np.array([np.bincount(c[c >= 0], minlength=4).argmax() for c in X.T]); cg = (cons[:-1] == 1) & (cons[1:] == 2)
        m = np.zeros(v3.L, bool); m[:-1] |= cg; m[1:] |= cg; Y = X.copy(); Y[:, m] = -1; return Y

    def oracle(X, dcols):
        Xd = X[:, dcols]; n = len(Xd); D = np.ones((n, n)); B = Xd >= 0
        for i in range(n):
            both = B[i] & B; m = both.sum(1); D[i] = np.where(m > 0, ((Xd[i] != Xd) & both).sum(1) / np.maximum(m, 1), 1.0)
        np.fill_diagonal(D, 0); return D

    res = []
    for seed in range(1, NSEED + 1):
        X, lab, dcols = simulate(seed); S = code_states(X); r = dict(seed=seed, diag_cols=len(dcols)); k = 6
        ev = lambda D: evaluate(D, lab, seed, k)[0]
        r["K2P"] = ev(v3.k2p_pd(X)); r["K2P CpG-masked"] = ev(v3.k2p_pd(cpg_mask(X)))
        ft = build_features(S); w, c = cov_weights(ft, S, seed); sel = w > 0
        selcols = set(int(x) for x in ft["col"][sel])
        r["selection"] = dict(features=int(ft["Fm"].shape[1]), selected=int(sel.sum()), threshold_z2=c,
                              selected_cols=len(selcols), selected_cols_true=len(selcols & set(dcols)),
                              true_cols_found=len(selcols & set(dcols)), true_cols=len(dcols))
        r["feat-all"] = ev(wham(ft["Fm"], np.ones(ft["Fm"].shape[1])))
        if sel.sum() >= 2:
            r["cov-selected"] = ev(wham(ft["Fm"][:, sel], np.ones(int(sel.sum()))))
            r["cov-weighted"] = ev(wham(ft["Fm"][:, sel], np.sqrt(w[sel])))
            for sh in (0.5, 0.75, 0.9): r[f"heavy{int(sh * 100)}"] = ev(wham(ft["Fm"], heavy(sel, sh)))
        fi = build_features(S, indel=True); wi, ci = cov_weights(fi, S, seed); si = wi > 0
        r["selection+indel"] = dict(features=int(fi["Fm"].shape[1]), selected=int(si.sum()), selected_indel=int((si & (fi["kind"] != "sub")).sum()))
        if si.sum() >= 2: r["cov-sel+indel"] = ev(wham(fi["Fm"][:, si], np.ones(int(si.sum()))))
        Z, kk = pca_pa(ft, seed); r["PCA-pa axes"] = kk
        if kk: r["PCA-pa"] = ev(euclid(Z))
        if dcols: r["oracle"] = ev(oracle(X, dcols))
        log = []; lv = nested(S, seed, log=log); nl = leaves_to_labels(lv, len(S))
        r["nested"] = dict(leaves=len(lv), ari=float(adjusted_rand_score(lab, nl)), sizes=sorted((int(len(v)) for v in lv.values()), reverse=True), log=log)
        res.append(r); json.dump(res, open(OUT, "w"), indent=1)
        print("seed", seed, json.dumps({a: b for a, b in r.items() if a != "nested"}), "\nnested", r["nested"]["leaves"], r["nested"]["ari"], r["nested"]["sizes"], flush=True)
    arms = [a for a in ["K2P", "K2P CpG-masked", "feat-all", "cov-selected", "cov-weighted", "heavy50", "heavy75", "heavy90", "cov-sel+indel", "PCA-pa", "oracle"] if any(a in r for r in res)]
    print(f"\n== DIAG={DIAG} NOSUB={NOSUB} TRUNC_FRAC={v3.TRUNC_FRAC}  mean ± sd over {len(res)} seeds (arm missing in a seed = not counted) ==")
    cols = ["knn_on_D", "silh_labels_on_D", "ari_umap_ward", "ari_tsne_ward", "hdb_ari"]
    print("%-16s" % "input" + "".join("%20s" % c for c in cols) + "   hdb_n per seed")
    for a in arms:
        v = [r[a] for r in res if a in r]
        print("%-16s" % a + "".join("%13.3f ±%.3f" % (np.mean([x[c] for x in v]), np.std([x[c] for x in v])) for c in cols) + "   " + str([x["hdb_n"] for x in v]) + f"  (n={len(v)})")
    print("selection:", [(r["selection"]["selected"], r["selection"]["selected_cols_true"], r["selection"]["selected_cols"], r["selection"]["true_cols"]) for r in res], "(selected features, selected cols that are true, selected cols, true cols)")
    print("selection+indel:", [(r["selection+indel"]["selected"], r["selection+indel"]["selected_indel"]) for r in res])
    print("PCA-pa axes:", [r["PCA-pa axes"] for r in res])
    print("nested: leaves", [r["nested"]["leaves"] for r in res], "ARI %.3f ±%.3f" % (np.mean([r["nested"]["ari"] for r in res]), np.std([r["nested"]["ari"] for r in res])))
