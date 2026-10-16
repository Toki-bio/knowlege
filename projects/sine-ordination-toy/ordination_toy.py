"""Toy: SINEderella-style copies with subfamily structure. PCA vs PCoA(K2P) vs NMDS vs t-SNE vs UMAP vs PaCMAP."""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import MDS, TSNE, trustworthiness
from sklearn.metrics import silhouette_score, adjusted_rand_score
from sklearn.cluster import AgglomerativeClustering
import umap, pacmap
rng=np.random.default_rng(1)
L=260; sub_names=list("ABCDEF")
# subfamily tree: root -> A ; A->B,C ; B->D,E ; C->F  (diagnostic mutations accumulate)
parent={"A":None,"B":"A","C":"A","D":"B","E":"B","F":"C"}
root=rng.integers(0,4,L); cons={}
for s in sub_names:
    base=root if parent[s] is None else cons[parent[s]].copy()
    n=0 if parent[s] is None else 6          # only ~6 diagnostic sites per step: subtle
    c=base.copy(); idx=rng.choice(L,n,replace=False)
    c[idx]=(c[idx]+rng.integers(1,4,n))%4; cons[s]=c
sizes=dict(A=250,B=300,C=250,D=250,E=200,F=150)
X=[];y=[]
for s in sub_names:
    for _ in range(sizes[s]):
        age=rng.gamma(2,0.045)                  # copy-level divergence (CpG/neutral decay)
        c=cons[s].copy(); m=rng.random(L)<age
        c[m]=(c[m]+rng.integers(1,4,m.sum()))%4
        X.append(c); y.append(s)
X=np.array(X); y=np.array(y); N=len(y)
lab=np.array([sub_names.index(v) for v in y])
# K2P distance (transitions: 0<->2, 1<->3 as A/G, C/T)
def k2p(X):
    D=np.zeros((len(X),len(X)))
    for i in range(len(X)):
        d=X[i]!=X; ts=((X[i]^2)==X)  # transition partner = x xor 2
        P=ts.mean(1); Q=(d&~ts).mean(1)
        a=1-2*P-Q; b=1-2*Q
        a=np.clip(a,1e-6,None); b=np.clip(b,1e-6,None)
        D[i]=-0.5*np.log(a)-0.25*np.log(b)
    np.fill_diagonal(D,0); return D
D=k2p(X)
OH=np.eye(4)[X].reshape(N,-1)
# 1 PCA
pca=PCA(10).fit(OH); Z_pca=pca.transform(OH)[:,:2]; ev_pca=pca.explained_variance_ratio_
# 2 PCoA on K2P (classical MDS)
J=np.eye(N)-1/N; B=-0.5*J@(D**2)@J; w,V=np.linalg.eigh(B); o=np.argsort(w)[::-1]; w,V=w[o],V[:,o]
pos=w[w>0]; ev_pcoa=w/pos.sum()    # share of positive-eigenvalue mass
Z_pcoa=V[:,:2]*np.sqrt(w[:2]); neg=-w[w<0].sum()/np.abs(w).sum()
# 3 NMDS (metric=False, init from PCoA)
Z_nmds=MDS(2,dissimilarity="precomputed",metric=False,n_init=1,init=None,random_state=1,max_iter=300,normalized_stress="auto").fit(D) if False else None
nm=MDS(2,metric=False,dissimilarity="precomputed",random_state=1,n_init=1,max_iter=300,init="random",normalized_stress="auto")
Z_nmds=nm.fit_transform(D); stress=nm.stress_
# 4-6 on PCA-30 of one-hot (standard pipeline) and t-SNE on K2P
P30=PCA(30,random_state=1).fit_transform(OH)
Z_tsne=TSNE(2,metric="precomputed",init="random",perplexity=40,random_state=1).fit_transform(D)
Z_umap=umap.UMAP(n_neighbors=30,min_dist=0.3,metric="precomputed",random_state=1).fit_transform(D)
Z_pacmap=pacmap.PaCMAP(n_components=2,n_neighbors=30,random_state=1).fit_transform(P30)
embs={"PCA (one-hot)":Z_pca,"PCoA (K2P)":Z_pcoa,"NMDS (K2P)":Z_nmds,"t-SNE (K2P)":Z_tsne,"UMAP (K2P)":Z_umap,"PaCMAP (PCA30)":Z_pacmap}
# metrics vs ground truth/original distances
def knn_label_acc(Z,k=15):
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.model_selection import cross_val_score
    return cross_val_score(KNeighborsClassifier(k),Z,lab,cv=5).mean()
rows=[]
for n,Z in embs.items():
    tw=trustworthiness(D,Z,n_neighbors=15,metric="precomputed") if False else trustworthiness(D,Z,n_neighbors=15)  # D as features is fine proxy? use precomputed below
    rows.append((n,silhouette_score(Z,lab),knn_label_acc(Z),adjusted_rand_score(lab,AgglomerativeClustering(6,linkage="ward").fit_predict(Z))))
# proper trustworthiness w/ precomputed original distances
from sklearn.manifold import trustworthiness as tw_
rows2=[]
for (n,s,k,a) in rows:
    t=tw_(D,embs[n],n_neighbors=15,metric="precomputed"); rows2.append((n,t,s,k,a))
print("%-16s %6s %6s %6s %6s"%("method","trust","silh","kNNacc","ARI"))
for r in rows2: print("%-16s %.3f %.3f %.3f %.3f"%r)
print("PCA PC1+PC2 var: %.1f%% (PC1 %.1f, PC2 %.1f); top10 %.1f%%"%(100*ev_pca[:2].sum(),100*ev_pca[0],100*ev_pca[1],100*ev_pca.sum()))
print("PCoA PC1+PC2 (of +eig mass): %.1f%%; negative-eigenvalue share %.1f%%"%(100*ev_pcoa[:2].sum(),100*neg))
print("NMDS stress (normalized): %.3f"%stress)
import json; json.dump(dict(rows=rows2,pca12=float(ev_pca[:2].sum()),pcoa12=float(ev_pcoa[:2].sum()),neg=float(neg),stress=float(stress)),open("ordination_metrics.json","w"))
# figure
cols=plt.cm.tab10(np.arange(6))
fig,ax=plt.subplots(2,4,figsize=(19,9.2))
for a,(n,Z) in zip(ax.flat[:6],embs.items()):
    for i,s in enumerate(sub_names): a.scatter(*Z[lab==i].T,s=6,color=cols[i],label=s,alpha=.7,lw=0)
    r=[r for r in rows2 if r[0]==n][0]
    ttl=n
    if n.startswith("PCA"): ttl+=f"\nPC1+PC2 = {100*ev_pca[:2].sum():.1f}% var"
    elif n.startswith("PCoA"): ttl+=f"\nAxes1+2 = {100*ev_pcoa[:2].sum():.1f}% of +eig mass"
    elif n.startswith("NMDS"): ttl+=f"\nstress = {stress:.3f}"
    else: ttl+="\naxes arbitrary (no % var)"
    a.set_title(f"{ttl}\ntrust={r[1]:.2f} kNN-acc={r[3]:.2f} ARI={r[4]:.2f}",fontsize=9.5); a.set_xticks([]);a.set_yticks([])
ax[0,0].legend(markerscale=3,title="subfamily",fontsize=8)
a=ax[1,2]; k=20
a.bar(np.arange(k)-.2,100*PCA(k).fit(OH).explained_variance_ratio_,.4,label="PCA (one-hot)")
a.bar(np.arange(k)+.2,100*ev_pcoa[:k],.4,label="PCoA (K2P)"); a.set_xlabel("axis"); a.set_ylabel("% variance"); a.set_title("Scree: why 2D PCA 'captures too little'\n(structure is spread over many weak axes)",fontsize=9.5); a.legend()
a=ax[1,3]; a.axis("off")
a.text(0,1,"Reading guide\n\n• PCA/PCoA axes: report % (PCoA: mind\n  negative eigenvalues; Cailliez/Lingoes\n  correction if large).\n• t-SNE/UMAP/PaCMAP: no % variance;\n  report trustworthiness, kNN-acc,\n  silhouette instead.\n• Don't read distances BETWEEN clusters\n  in t-SNE/UMAP. Sweep n_neighbors /\n  perplexity; trust only stable clusters.\n• Subfamilies here differ by ~6 diag.\n  sites on top of ~9% mean copy noise.",va="top",fontsize=10,family="monospace")
fig.suptitle("Toy SINEderella output: 1,400 copies, 6 subfamilies (tree-structured), 260 bp — ordination comparison",fontsize=13)
plt.tight_layout(); plt.savefig("fig1_ordination_comparison.png",dpi=150)
# sweep stability for UMAP
fig,ax=plt.subplots(1,4,figsize=(16,4))
for a,nn in zip(ax,[5,15,50,150]):
    Z=umap.UMAP(n_neighbors=nn,min_dist=.3,metric="precomputed",random_state=1).fit_transform(D)
    for i,s in enumerate(sub_names): a.scatter(*Z[lab==i].T,s=5,color=cols[i],lw=0)
    a.set_title(f"UMAP n_neighbors={nn}\ntrust={tw_(D,Z,n_neighbors=15,metric='precomputed'):.2f}"); a.set_xticks([]);a.set_yticks([])
plt.suptitle("Parameter sweep (always do this): stable structure only"); plt.tight_layout(); plt.savefig("fig2_umap_sweep.png",dpi=140)
