---
title: Beyond PCA for SINEderella output; pan-SINEome visuals
type: notes
---

> Type: notes

# Beyond PCA for SINEderella output, and visuals for orthologous loci / pan-SINEome

Toy-data study (simulated, **not** real genomes). Code, simulated matrix and metrics: [`projects/sine-ordination-toy/`](https://github.com/Toki-bio/knowlege/tree/main/projects/sine-ordination-toy). Figures are re-rendered by running the two scripts (Python 3.13, numpy, scipy, scikit-learn, matplotlib, umap-learn 0.5.12, pacmap 0.9.1).

Related: [dimensionality reduction for SINE separation](../dim-reduction-for-sine-separation/), [clustering for SINE subfamily discovery](../clustering-for-sine-subfamily-discovery/), [TE-anchor graph pan-SINEome](../../projects/te-anchor-graph-pansineome/), [sq2 distant orthologous loci](../../projects/sq2-distant-orthologous-loci/).

## 1. Problem

PCA of SINEderella copies explains little variance in PC1+PC2, so the 2D plot looks muddy and the percentage on the axes looks bad.

## 2. Toy ordination experiment

Simulation (`ordination_toy.py`, seed fixed): 1,400 copies, 260 bp ungapped, 6 subfamilies on a tree (A→B,C; B→D,E; C→F), only ~6 diagnostic substitutions per branch, plus per-copy noise (gamma-distributed divergence, mean ≈ 9%). Features: one-hot (1,040 columns) or K2P distance matrix. Metrics are computed against the true subfamily labels and the K2P distances.

| Method | trustworthiness (k=15) | silhouette | kNN acc | ARI (Ward, k=6) |
|---|---|---|---|---|
| PCA (one-hot) | 0.687 | 0.552 | 0.965 | 0.905 |
| PCoA (K2P) | 0.696 | 0.357 | 0.886 | 0.592 |
| NMDS (K2P) | 0.681 | 0.052 | 0.901 | 0.283 |
| t-SNE (K2P) | 0.699 | 0.755 | 0.996 | 0.993 |
| UMAP (K2P) | 0.679 | 0.860 | 0.997 | 0.993 |
| PaCMAP (PCA-30) | 0.678 | 0.893 | 0.992 | 0.981 |

Variance / stress: PCA PC1+PC2 = **11.9 %** (PC1 8.7, PC2 3.2; top 10 PCs 21.1 %). PCoA axes 1+2 = 12.1 % of the positive-eigenvalue mass, with a **33 % negative-eigenvalue share** (the K2P matrix is far from Euclidean at this noise level). NMDS normalized stress 0.251.

![Ordination comparison]({{ '/assets/img/sine-ordination/fig1_ordination_comparison.png' | relative_url }})

![UMAP n_neighbors sweep]({{ '/assets/img/sine-ordination/fig2_umap_sweep.png' | relative_url }})

### Reading

- The low variance is not a bug: structure is spread over many weak axes (scree in the figure), so the *percentage* is the wrong yardstick, not necessarily PCA's picture.
- **PCoA on K2P and NMDS did not beat PCA here.** PCoA only helps when you need a % explained for a metric distance, and then needs a negative-eigenvalue correction.
- **t-SNE, UMAP and PaCMAP recovered the subfamilies** (ARI 0.98–0.99). Trustworthiness is ~0.68–0.70 for every method, so it did not discriminate in this noisy toy.
- UMAP/t-SNE/PaCMAP give no % variance; report kNN accuracy / silhouette / trustworthiness, always show a parameter sweep, and do not interpret distances between clusters.
- Practical pipeline to try on real data: K2P matrix (or PCA-30 of one-hot) → UMAP or PaCMAP → cluster (HDBSCAN/Leiden), keep the PCA/PCoA scree as a side panel.
- Limits: no indels, one noise model, one seed; hence only a demonstration. Needs to be repeated on real SINEderella output before any conclusion.

## 3. Toy pan-SINEome and orthologous-locus visuals

Simulation (`pansineome_toy.py`): 36 species in 4 clades (Gekkota 8, Iguania 10, Serpentes 10, Lacertoidea 8), random tree, 3,500 loci inserted irreversibly on branches (some branches with 6× insertion bursts), 3 % of present cells excised, 1.5 % marked missing / flank-ambiguous. 3,444 loci retained. Matrix: `toy_pansineome_PAV.tsv`.

**Locus × species heatmap on the species tree** (the main figure; three states: present / empty site / missing):

![Tree + presence/absence heatmap]({{ '/assets/img/sine-ordination/fig3_tree_heatmap.png' | relative_url }})

**Composition, per-species counts, accumulation curves with Heaps' law fit** (α = 0.40 on the toy):

![Pan-SINEome summary]({{ '/assets/img/sine-ordination/fig4_pansineome_summary.png' | relative_url }})

**Clade-level UpSet** (hand-written in matplotlib; `upsetplot` 0.9.0 failed on pandas 3):

![UpSet by clade]({{ '/assets/img/sine-ordination/fig4b_upset_clades.png' | relative_url }})

**Jaccard PCoA on presence/absence** (Jaccard is a metric, so eigenvalue % is meaningful; PCoA1+2 = 71.4 % on the toy) and **UPGMA recovering the clades**:

![Jaccard PCoA]({{ '/assets/img/sine-ordination/fig5_jaccard_pcoa.png' | relative_url }})

![UPGMA on PAV]({{ '/assets/img/sine-ordination/fig7_pav_upgma.png' | relative_url }})

**Orthologous-locus views** (one locus as flank/SINE/flank strip across species in tree order, locus-quality scatter, frequency spectrum). Raster-based, so it avoids per-character DOM spans (the known slowness of the browser SINE Loci Alignment Viewer):

![Ortholog locus views]({{ '/assets/img/sine-ordination/fig6_ortholog_locus_views.png' | relative_url }})

### Toy caveats

- Core = 0 here by construction (no insertions on the root branch); not a biological result.
- Flank-identity values in the locus-quality panel are simulated.
- Heaps' α depends on the phylogeny sampling; treat open/closed calls with care since species are not independent.
- Insertion presence/absence is near-homoplasy-free but not perfect (ILS, hotspots, excision), which is why the toy includes excision and missing cells.

## 4. Literature pointers (from the search round; links not re-verified)

- UMAP / t-SNE / PaCMAP vs PCA: axes arbitrary, inter-cluster distance not meaningful, can create artificial clusters; use for cluster discovery, sweep n_neighbors/perplexity.
- TE subfamily clustering precedents: DBSCAN on K2P from MSA (MCHelper), SCULU, k-mer classifiers (TEclass, DeepTE, NeuralTE, BERTE), similarity network + Leiden.
- Pangenome-style visuals: tree + heatmap, UpSet (upsetplot), core/shell/cloud/singleton, pan/core accumulation with Heaps' law, Panache, PanViz.
- SINE presence/absence as phylogenetic markers: Myotis (Ves SINE insertion matrix), SplitsTree networks, Camin–Sokal parsimony.

## 5. Next steps

1. Run the ordination script on a real SINEderella output (copies of one SINE family, one MSA) and compare PCA vs UMAP/PaCMAP with the same metrics.
2. Build the real sq2 presence/absence matrix (SINE_orth_loc output) and render the tree + heatmap and the clade UpSet.
3. Port the locus strip (fig 6) to the MSA/SINE viewer via canvas.
