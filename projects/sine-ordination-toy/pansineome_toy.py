"""Toy: pan-SINEome (orthologous-locus presence/absence) on a species tree + ortholog-locus views."""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from scipy.optimize import curve_fit
from sklearn.manifold import MDS
from upsetplot import UpSet, from_memberships
rng=np.random.default_rng(7)
# ---- random species tree (nested lists), 36 species in 4 clades, branch lengths
class Node:
    def __init__(s,name=None,kids=(),bl=0): s.name,s.kids,s.bl=name,list(kids),bl
def build(names,depth=0):
    if len(names)==1: return Node(names[0],bl=rng.uniform(.3,1.2))
    k=rng.integers(1,len(names)); 
    return Node(kids=[build(names[:k],depth+1),build(names[k:],depth+1)],bl=rng.uniform(.4,1.5)/(1+depth*.15))
clades={"Gekkota":8,"Iguania":10,"Serpentes":10,"Lacertoidea":8}
sp=[];root=Node(); 
for c,n in clades.items():
    names=[f"{c[:3]}_{i+1:02d}" for i in range(n)]; sp+=names; t=build(names); t.bl=3.0; root.kids.append(t)
# fold clades into a caterpillar backbone
def fold(k):
    if len(k)==1: return k[0]
    return Node(kids=[k[0],fold(k[1:])],bl=2.0)
root=fold(root.kids); root.bl=0
# ---- insertions: each locus inserted on one branch (rate ~ length*hotspot); irreversible; deletions/missing noise
nodes=[];
def walk(n,anc):
    n.id=len(nodes); nodes.append(n); n.anc=anc
    for k in n.kids: walk(k,n.id)
walk(root,-1)
def leaves(n): return [n.name] if not n.kids else sum((leaves(k) for k in n.kids),[])
leafset=[leaves(n) for n in nodes]; si={s:i for i,s in enumerate(sp)}
w=np.array([n.bl for n in nodes]); w[0]=0
w=w*np.where(rng.random(len(w))<.15,6,1)   # insertion bursts on some branches (subfamily bursts)
w/=w.sum(); LOCI=3500
br=rng.choice(len(nodes),LOCI,p=w)
M=np.zeros((LOCI,len(sp)),int)
for j,b in enumerate(br):
    M[j,[si[s] for s in leafset[b]]]=1
true_branch=br.copy()
# homoplasy-like: precise excision in 3% of present cells; and 1.5% "missing" (no usable flank)
exc=(M==1)&(rng.random(M.shape)<.03); M[exc]=0
miss=rng.random(M.shape)<.015
Mobs=M.copy(); status=np.where(miss,2,M)   # 0 empty site, 1 present, 2 missing/ambiguous
M=np.where(miss,0,M)
# drop loci present nowhere
keep=M.sum(1)>0; M=M[keep];status=status[keep];true_branch=true_branch[keep]
print("loci",M.shape)
# ---- leaf order from tree
order=[s for s in leaves(root)]; oi=[si[s] for s in order]
# ---- Fig3: tree + heatmap
def xy(n,depth=0):
    if not n.kids: n.y=oi.index(si[n.name]); n.x=depth+n.bl; return
    for k in n.kids: xy(k,depth+n.bl)
    n.y=np.mean([k.y for k in n.kids]); n.x=depth+n.bl
xy(root)
def draw(n,ax,px=0):
    ax.plot([px,n.x],[n.y,n.y],"k",lw=.8)
    if n.kids:
        ys=[k.y for k in n.kids]; ax.plot([n.x,n.x],[min(ys),max(ys)],"k",lw=.8)
        for k in n.kids: draw(k,ax,n.x)
# sort loci by (first present species in tree order, count) -> block-diagonal structure
Mo=M[:,oi]; So=status[:,oi]
key=sorted(range(len(Mo)),key=lambda j:(-Mo[j].sum() if Mo[j].sum()>=len(sp)*.9 else 0, tuple(-Mo[j])))
Mo2=Mo[key]; So2=So[key]; sums=Mo2.sum(1)
core=sums>=0.99*len(sp)*0+0.95*len(sp)  # >=95% species
shell=(sums>=0.15*len(sp))&~core; cloud=(sums>1)&~core&~shell; single=sums==1
fig=plt.figure(figsize=(14,9)); gs=fig.add_gridspec(1,3,width_ratios=[1.1,4,.01],wspace=.02)
a0=fig.add_subplot(gs[0]); draw(root,a0); a0.set_ylim(len(sp)-.5,-.5); a0.axis("off")
a1=fig.add_subplot(gs[1]); cm=ListedColormap(["#f2f2f2","#1b4f72","#e67e22"])
a1.imshow(So2.T,aspect="auto",cmap=cm,interpolation="nearest",vmin=0,vmax=2)
a1.set_yticks(range(len(sp))); a1.set_yticklabels(order,fontsize=7); a1.yaxis.tick_right(); a1.set_xlabel(f"{len(Mo2)} orthologous SINE loci (sorted by phylogenetic distribution)")
a1.legend(handles=[Patch(color="#1b4f72",label="SINE present"),Patch(color="#f2f2f2",ec="k",label="empty site (absent)"),Patch(color="#e67e22",label="missing / flank ambiguous")],loc="upper center",bbox_to_anchor=(.5,1.07),ncol=3,frameon=False)
fig.suptitle("Pan-SINEome: locus × species presence/absence aligned to species tree",y=.98); plt.savefig("fig3_tree_heatmap.png",dpi=150,bbox_inches="tight")
# ---- Fig4: composition + UpSet + accumulation
cat=np.select([single,core,shell,cloud],["singleton","core (≥95%)","shell (15–95%)","cloud (2 spp–15%)"],"other")
# clade-level UpSet
cl={c:[i for i,s in enumerate(sp) if s.startswith(c[:3])] for c in clades}
memb=[tuple(c for c in clades if M[j,cl[c]].mean()>=.5) for j in range(len(M))]
memb=[m for m in memb if m]
ups=from_memberships(memb); 
fig=plt.figure(figsize=(17,4.8)); gs=fig.add_gridspec(1,3,wspace=.3)
# composition bars by clade
ax=fig.add_subplot(gs[0]); cats=["core (≥95%)","shell (15–95%)","cloud (2 spp–15%)","singleton"]
v=[(cat==c).sum() for c in cats]; ax.barh(cats[::-1],v[::-1],color=["#c0392b","#e67e22","#2e86c1","#1b4f72"]); ax.set_title("Pan-SINEome composition (loci)")
for i,x in enumerate(v[::-1]): ax.text(x,i,f" {x}",va="center")
# species-wise private loci
ax=fig.add_subplot(gs[1]); priv=[(M[single][:,i]).sum() for i in range(len(sp))]; tot=M.sum(0)
ax.bar(range(len(sp)),tot[oi],color="#1b4f72",label="all loci"); ax.bar(range(len(sp)),[priv[i] for i in oi],color="#e74c3c",label="species-private")
ax.set_xticks(range(len(sp))); ax.set_xticklabels(order,rotation=90,fontsize=6); ax.legend(); ax.set_title("SINE loci per species")
# accumulation
def acc(M,n=60):
    P=[];C=[]
    for _ in range(n):
        p=rng.permutation(M.shape[1]); cum=np.maximum.accumulate(M[:,p],1)  # not right for core; do explicit
        seen=np.cumsum(M[:,p],1); P.append((seen>0).sum(0)); C.append((seen==np.arange(1,M.shape[1]+1)).sum(0))
    return np.array(P),np.array(C)
P,C=acc(M); x=np.arange(1,len(sp)+1)
heaps=lambda n,k,a:k*n**a
(k,al),_=curve_fit(heaps,x,P.mean(0),p0=[500,.3],maxfev=10000)
ax=fig.add_subplot(gs[2]); ax.fill_between(x,P.min(0),P.max(0),alpha=.2,color="#1b4f72"); ax.plot(x,P.mean(0),"o-",ms=3,color="#1b4f72",label="pan (all loci)")
ax.fill_between(x,C.min(0),C.max(0),alpha=.2,color="#c0392b"); ax.plot(x,C.mean(0),"s-",ms=3,color="#c0392b",label="core (all genomes)")
ax.plot(x,heaps(x,k,al),"k--",lw=1,label=f"Heaps' law fit  n^α, α={al:.2f}"); ax.set_xlabel("genomes added (60 random orders)"); ax.set_ylabel("loci"); ax.legend(fontsize=8)
ax.set_title("Accumulation: α<1 ≈ open pan-SINEome\n(α is OPEN-limit caveat: phylogeny-structured)",fontsize=9)
plt.savefig("fig4_pansineome_summary.png",dpi=140,bbox_inches="tight")
from collections import Counter
cnt=Counter(memb).most_common(14); names=list(clades)
fig=plt.figure(figsize=(10,5.5)); g=fig.add_gridspec(2,2,width_ratios=[1,3],height_ratios=[2.2,1.3],hspace=.05,wspace=.05)
ab=fig.add_subplot(g[0,1]); ab.bar(range(len(cnt)),[c for _,c in cnt],color="#1b4f72"); ab.set_xticks([])
for i,(_,c) in enumerate(cnt): ab.text(i,c,str(c),ha="center",va="bottom",fontsize=8)
ab.set_ylabel("loci"); ab.set_title("Clade-level UpSet (own matplotlib impl.; locus assigned to a clade if >=50% of its species carry it)",fontsize=9)
am=fig.add_subplot(g[1,1],sharex=ab)
for i,(m,_) in enumerate(cnt):
    ys=[names.index(c) for c in m]; am.plot([i]*len(ys),ys,"-",color="k",lw=2)
    for r in range(len(names)): am.scatter(i,r,s=60,color="k" if names[r] in m else "#dddddd",zorder=3)
am.set_yticks(range(len(names))); am.set_yticklabels(names); am.set_ylim(-.5,len(names)-.5); am.invert_yaxis()
axl=fig.add_subplot(g[1,0]); axl.barh(range(len(names)),[sum(names[r] in m for m in memb) for r in range(len(names))],color="#888"); axl.invert_xaxis(); axl.invert_yaxis(); axl.set_yticks([]); axl.set_xlabel("loci per clade")
plt.savefig("fig4b_upset_clades.png",dpi=140,bbox_inches="tight")
print("categories",dict(zip(cats,v)),"heaps alpha %.3f k %.1f"%(al,k))
# ---- Fig5: Jaccard PCoA + distance heatmap
inter=M.T@M; un=M.sum(0)[:,None]+M.sum(0)[None]-inter; Dj=1-inter/un
Jm=np.eye(len(sp))-1/len(sp); B=-.5*Jm@Dj**2@Jm; ww,V=np.linalg.eigh(B); o=np.argsort(ww)[::-1]; ww,V=ww[o],V[:,o]
evp=ww/ww[ww>0].sum(); Z=V[:,:2]*np.sqrt(ww[:2])
fig,ax=plt.subplots(1,2,figsize=(13,5.5)); cc=dict(zip(clades,["#1b9e77","#d95f02","#7570b3","#e7298a"]))
for i,s in enumerate(sp): ax[0].scatter(*Z[i],color=cc[[c for c in clades if s.startswith(c[:3])][0]],s=40)
for c in clades: ax[0].scatter([],[],color=cc[c],label=c)
ax[0].legend(); ax[0].set_xlabel(f"PCoA1 ({100*evp[0]:.1f}%)"); ax[0].set_ylabel(f"PCoA2 ({100*evp[1]:.1f}%)"); ax[0].set_title("PCoA of Jaccard distance on SINE presence/absence\n(metric ⇒ eigenvalue % meaningful)")
im=ax[1].imshow(Dj[np.ix_(oi,oi)],cmap="viridis_r"); ax[1].set_title("Jaccard distance (tree order)"); plt.colorbar(im,ax=ax[1],fraction=.046)
plt.tight_layout(); plt.savefig("fig5_jaccard_pcoa.png",dpi=150)
print("PCoA1+2 %.1f%%"%(100*evp[:2].sum()))
np.savetxt("toy_pansineome_PAV.tsv",M,fmt="%d",delimiter="\t",header="\t".join(sp),comments="")
# ---- Fig6: ortholog-locus panels (sq2-like): per-locus alignment strip + locus QC view
nL=14; picks=[j for j in rng.choice(np.where(shell|cloud)[0],nL,replace=False)]
fl=40; sl=60; Lw=2*fl+sl
fig,axs=plt.subplots(1,3,figsize=(17,7),gridspec_kw=dict(width_ratios=[3,1.3,1.3]))
ax=axs[0]; j=picks[0]; img=np.zeros((len(sp),Lw)); 
for r,i in enumerate(oi):
    flank_id=np.clip(1-.04*abs(rng.normal(0,1))*(1+r/6),0,1)   # flank conservation decays with distance
    img[r,:fl]=rng.random(fl)<flank_id; img[r,fl+sl:]=rng.random(fl)<flank_id
    if So2[0,r]==2: img[r]=.5
    elif Mo2[0,r]==1 or True:
        pass
jj=key.index(np.where(~single)[0][0]) if False else None
# choose a mid-frequency locus
cand=np.where((sums>8)&(sums<24))[0]; ji=cand[len(cand)//2]
for r in range(len(sp)):
    i=oi[r]; fid=np.clip(.995-.0045*r,0,1)
    row=np.zeros(Lw); row[:fl]=(rng.random(fl)<fid)*1; row[fl+sl:]=(rng.random(fl)<fid)*1
    if So2[ji,r]==1: row[fl:fl+sl]=2+(rng.random(sl)<.9)*0.5   # SINE body
    elif So2[ji,r]==0: row[fl:fl+sl]=-1                          # empty site: flanks joined (gap in alignment)
    else: row[:]=-2
    img[r]=row
cm2=ListedColormap(["#bbbbbb","#ffffff","#e8b4b4","#1b4f72","#2e86c1"])  # -2 missing,-1 gap,0 mismatch,1 match,2..SINE
import matplotlib.colors as mc
disp=np.select([img==-2,img==-1,img==0,img==1,img>=2],[0,1,2,3,4],0)
ax.imshow(disp,aspect="auto",cmap=cm2,interpolation="nearest",vmin=0,vmax=4)
ax.axvline(fl-.5,color="k",lw=.6);ax.axvline(fl+sl-.5,color="k",lw=.6)
ax.set_yticks(range(len(sp)));ax.set_yticklabels(order,fontsize=6);ax.set_xticks([fl/2,fl+sl/2,fl+sl+fl/2]);ax.set_xticklabels(["5' flank","SINE","3' flank"])
ax.legend(handles=[Patch(color="#1b4f72",label="flank match"),Patch(color="#e8b4b4",label="flank mismatch"),Patch(color="#2e86c1",label="SINE insert"),Patch(color="#ffffff",ec="k",label="empty site (gap)"),Patch(color="#bbbbbb",label="missing")],fontsize=7,loc="lower center",bbox_to_anchor=(.5,-.17),ncol=5,frameon=False)
ax.set_title(f"One orthologous locus across {len(sp)} genomes (aligned to tree)\nraster/SVG, not per-character DOM")
# locus QC: flank identity (mean of 5'/3') vs # species carrying, color by category
ax=axs[1]; fi=np.clip(rng.normal(.93,.05,len(M)),.5,1)-.0*sums
ax.scatter(sums,fi,s=3,alpha=.4,c=np.select([single,core,shell,cloud],["#1b4f72","#c0392b","#e67e22","#2e86c1"],"#999999"))
ax.axhline(.85,color="k",ls="--",lw=.8); ax.set_xlabel("species carrying SINE"); ax.set_ylabel("min flank identity"); ax.set_title("Locus QC (toy values)\nfilter < threshold → 'ambiguous'",fontsize=9)
# insertion-age proxy: carrying-clade-size, distribution of locus "tree depth" (clade-specific count)
ax=axs[2]; ax.hist(sums,bins=np.arange(1,len(sp)+2)-.5,color="#1b4f72"); ax.set_yscale("log"); ax.set_xlabel("# species carrying locus"); ax.set_ylabel("loci (log)"); ax.set_title("Frequency spectrum (SFS-like)\nU-shape ⇒ many singletons + core",fontsize=9)
plt.tight_layout(); plt.savefig("fig6_ortholog_locus_views.png",dpi=150,bbox_inches="tight")
# ---- Fig7: SINE insertion phylogeny via parsimony-informative sites -> NJ-like check: compare Jaccard UPGMA to true tree order
from scipy.cluster.hierarchy import linkage,dendrogram
from scipy.spatial.distance import squareform
Zl=linkage(squareform(Dj,checks=False),"average"); fig,ax=plt.subplots(figsize=(11,4)); dendrogram(Zl,labels=sp,ax=ax,leaf_font_size=7); ax.set_title("UPGMA on Jaccard(SINE PAV) recovers clades"); plt.tight_layout(); plt.savefig("fig7_pav_upgma.png",dpi=140)
