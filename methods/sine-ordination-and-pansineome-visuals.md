---
title: Beyond PCA for SINEderella output; pan-SINEome visuals
type: notes
---

> Type: notes

# Beyond PCA for SINEderella output, and visuals for orthologous loci / pan-SINEome

**Status: toy data only.** Everything on this page comes from two simulations, not from real genomes. It shows what the plots look like and which numbers to report. It supports no conclusion about any real SINE family or any real set of loci.

Code, simulated matrix and metrics: [`projects/sine-ordination-toy/`](https://github.com/Toki-bio/knowlege/tree/main/projects/sine-ordination-toy).

Related: [dimensionality reduction for SINE separation](../dim-reduction-for-sine-separation/), [clustering for SINE subfamily discovery](../clustering-for-sine-subfamily-discovery/), [HDBSCAN vs vsearch for SINE grouping](../hdbscan-vs-vsearch-for-sine-grouping/), [TE-insertion anchor graph](../../projects/te-anchor-graph-pansineome/), [sq2 distant orthologous loci](../../projects/sq2-distant-orthologous-loci/).

## In short

- A 2D PCA of SINE copies shows a low percentage on its axes because the subfamily signal is spread over many weak axes. The low percentage is expected, and by itself it does not make the PCA picture wrong.
- On the toy, t-SNE, UMAP and PaCMAP drew the six simulated subfamilies as separate groups. PCA drew them as overlapping groups. PCoA and NMDS on K2P distances scored below PCA on silhouette, kNN accuracy and ARI.
- The comparison is weaker than the table suggests: three of the four scores are measured inside the 2D picture, the methods did not all get the same input, NMDS was run with one random start, and the toy has no indels and one seed. See [2.5](#25-what-the-toy-cannot-show).
- For presence/absence of orthologous loci, the tree + heatmap with three states (present, empty site, missing) is the main figure. The other panels are summaries of it.
- No real presence/absence matrix exists yet, so every real-data item in [section 4](#4-further-developments) is a proposal.

## 1. The problem

PCA of SINEderella copies explains little variance in PC1 + PC2. The 2D plot looks muddy and the percentages on the axes look bad. The question was whether another ordination gives a clearer picture, and what to report instead of a variance percentage.

A reminder from the related pages: an ordination is a picture. It does not assign copies to subfamilies.

## 2. Toy ordination experiment

Script: `ordination_toy.py`.

### 2.1 What was simulated

| Item | Value |
|---|---|
| Copies | 1,400, each 260 bp, no gaps |
| Subfamilies | 6 (A 250, B 300, C 250, D 250, E 200, F 150 copies) |
| Subfamily tree | A is the root; B and C derive from A; D and E from B; F from C |
| Diagnostic sites | 6 substitutions at each step of the tree |
| Copy-level noise | each copy has its own substitution rate, drawn from a gamma distribution (shape 2, scale 0.045; mean 9 % of sites) |
| Random seed | 1 |
| Features | one-hot matrix (1,040 columns), or K2P distance matrix |

### 2.2 Methods and settings

| Method | Input | Settings in the script |
|---|---|---|
| PCA | one-hot | 10 components, first two plotted |
| PCoA | K2P | classical scaling, no negative-eigenvalue correction |
| NMDS | K2P | scikit-learn `MDS(metric=False)`, one random start, 300 iterations |
| t-SNE | K2P | perplexity 40, random initialization |
| UMAP | K2P | `n_neighbors` 30, `min_dist` 0.3, default initialization |
| PaCMAP | first 30 PCs of one-hot | `n_neighbors` 30, other settings default |

The UMAP sweep figure uses `n_neighbors` 5, 15, 50 and 150.

### 2.3 Results

| Method | Trustworthiness | Silhouette | kNN accuracy | ARI |
|---|---|---|---|---|
| PCA (one-hot) | 0.687 | 0.552 | 0.965 | 0.905 |
| PCoA (K2P) | 0.696 | 0.357 | 0.886 | 0.592 |
| NMDS (K2P) | 0.681 | 0.052 | 0.901 | 0.283 |
| t-SNE (K2P) | 0.699 | 0.755 | 0.996 | 0.993 |
| UMAP (K2P) | 0.679 | 0.860 | 0.997 | 0.993 |
| PaCMAP (PCA-30) | 0.678 | 0.893 | 0.992 | 0.981 |

Where each score is measured:

| Score | Measured on | Meaning |
|---|---|---|
| Trustworthiness (k = 15) | K2P distances against the 2D picture | are the 15 nearest neighbors in the picture also close in sequence distance |
| Silhouette | the 2D picture, with the true subfamily labels | how compact and separate the true subfamilies look in the picture |
| kNN accuracy (k = 15, 5-fold) | the 2D picture, with the true labels | can a copy's subfamily be read from its neighbors in the picture |
| ARI | Ward clustering (k = 6) of the 2D picture, against the true labels | does clustering the picture return the true subfamilies |

Variance and stress:

- PCA: PC1 + PC2 = **11.9 %** (PC1 8.7, PC2 3.2); the first 10 PCs together 21.1 %.
- PCoA: axes 1 + 2 = 12.1 % of the positive-eigenvalue mass. Negative eigenvalues make up **33 %** of the summed absolute eigenvalues, so the K2P matrix is far from Euclidean at this noise level.
- NMDS: normalized stress 0.251.

![Ordination comparison]({{ '/assets/img/sine-ordination/fig1_ordination_comparison.png' | relative_url }})

![UMAP n_neighbors sweep]({{ '/assets/img/sine-ordination/fig2_umap_sweep.png' | relative_url }})

### 2.4 What the toy shows

- **The low variance is not a bug.** The scree panel shows the structure spread over many weak axes. The percentage is the wrong yardstick for a 2D plot of such data.
- **PCA already separates most of it** (ARI 0.905), with the neighboring subfamilies overlapping.
- **t-SNE, UMAP and PaCMAP separate all six subfamilies** (ARI 0.981 to 0.993).
- **Trustworthiness is 0.68 to 0.70 for every method**, so it did not tell the methods apart here.
- UMAP, t-SNE and PaCMAP have no variance percentage, and their axes have no meaning. The position of one group relative to another should not be read.

### 2.5 What the toy cannot show

- **Three of the four scores are measured in the 2D picture.** Silhouette, kNN accuracy and ARI all reward a method for drawing compact, separate groups. That is what UMAP and PaCMAP are built to do, so their high silhouette (0.86, 0.89) is partly a property of the method. Only trustworthiness refers back to sequence distances, and it does not separate the methods.
- **The methods did not get the same input.** PCA used one-hot, PaCMAP used the first 30 PCs, the others used K2P. Method and input cannot be told apart in this table.
- **NMDS was run weakly.** One random start and 300 iterations; the NMDS panel is star-shaped, which looks like a run that did not converge. "NMDS did not beat PCA" may describe this run and not the method.
- **Initialization differs.** t-SNE started from random positions and UMAP from its default. The large-scale arrangement of groups depends on initialization (Kobak & Linderman 2021), so the arrangements in the figure are not comparable.
- **The toy is easy.** No indels, no truncated copies, one noise model, one seed. An ARI of 0.99 leaves no room to compare anything better.
- **K is given.** ARI uses Ward clustering with the true number of subfamilies (6). On real data that number is unknown.

## 3. Toy pan-SINEome and orthologous-locus visuals

Script: `pansineome_toy.py`. Matrix: `toy_pansineome_PAV.tsv`.

### 3.1 What was simulated

| Item | Value |
|---|---|
| Species | 36 in 4 groups named Gekkota (8), Iguania (10), Serpentes (10), Lacertoidea (8) |
| Tree | random inside each group; groups joined as a ladder. It is a toy tree, not a squamate phylogeny |
| Loci | 3,500, each inserted once on one branch and inherited by all descendants |
| Insertion rate | proportional to branch length; 6 times higher on some branches (each branch with probability 0.15); none on the root branch |
| Excision | 3 % of "present" cells set to absent |
| Missing | 1.5 % of all cells marked missing / flank-ambiguous |
| Loci kept | 3,444 (present in at least one species) |
| Random seed | 7 |

Missing cells are drawn as a separate state in the heatmap and in the locus strip. In every count and distance below they are treated as absent.

### 3.2 The figures

**Tree + heatmap.** One column per locus, one row per species in tree order, three states: present, empty site, missing. This is the main figure; the others summarize it.

![Tree + presence/absence heatmap]({{ '/assets/img/sine-ordination/fig3_tree_heatmap.png' | relative_url }})

**Composition, loci per species, accumulation.**

| Category | Definition | Loci |
|---|---|---|
| Core | in at least 95 % of species | 0 |
| Shell | in 15 % to 95 % of species | 1,266 |
| Cloud | in 2 species up to 15 % | 1,001 |
| Singleton | in 1 species | 1,177 |

The accumulation curve is the mean over 60 random orders of adding genomes. The fit of k·n^α to it gives α = 0.40.

![Pan-SINEome summary]({{ '/assets/img/sine-ordination/fig4_pansineome_summary.png' | relative_url }})

**UpSet by group.** A locus counts for a group when at least half of that group's species carry it. The 1,863 loci that reach this in no group (all singletons, among others) are left out. Only 6 combinations occur, and all are drawn. The plot is hand-written in matplotlib because `upsetplot` 0.9.0 failed on pandas 3.

![UpSet by clade]({{ '/assets/img/sine-ordination/fig4b_upset_clades.png' | relative_url }})

**Jaccard PCoA.** PCoA 1 + 2 = 71.4 %. The script does not print the negative-eigenvalue share; recomputed from `toy_pansineome_PAV.tsv`, this matrix has no negative eigenvalues, so the percentage is of the whole.

![Jaccard PCoA]({{ '/assets/img/sine-ordination/fig5_jaccard_pcoa.png' | relative_url }})

**UPGMA on Jaccard distance.** The dendrogram shows the four groups and the same order of joining as the simulated tree. This was judged by eye; no tree distance was computed.

![UPGMA on PAV]({{ '/assets/img/sine-ordination/fig7_pav_upgma.png' | relative_url }})

**Orthologous-locus views.** Left: one locus as a 5′ flank / SINE / 3′ flank strip across all species in tree order. It is drawn as a raster, which avoids one DOM element per character (the known slowness of the browser SINE Loci Alignment Viewer). Middle: locus quality scatter. Right: number of species per locus.

![Ortholog locus views]({{ '/assets/img/sine-ordination/fig6_ortholog_locus_views.png' | relative_url }})

### 3.3 Toy caveats

- **Core = 0 by construction.** No insertions were placed on the root branch. It is not a biological result.
- **Flank values are invented.** In the strip, flank matches are random draws with a match probability that simply falls from the top row to the bottom row. In the quality scatter, flank identities are random numbers not linked to the loci.
- **The spectrum panel title is wrong for this toy.** It says "U-shape ⇒ many singletons + core", but there is no core. The peaks are the sizes of the simulated groups.
- **Heaps' α is descriptive only.** It depends on which species are sampled, and species on a tree are not independent.
- **PCoA percentages need the negative-eigenvalue share next to them.** Here the Jaccard matrix has none, but that was checked afterwards, not known in advance. On a real matrix the share should be printed with the percentage.
- **Presence/absence is near homoplasy-free, not homoplasy-free** (incomplete lineage sorting, hot spots, excision). This is why the toy has excision and missing cells.

## 4. Further developments

Rules that every item below has to respect:

- Only copies with tested-unique flanks are used.
- Flanks are compared as short pieces. No step aligns long stretches or whole genomes.
- A target-site duplication is a bonus fact, never a requirement. Old insertions often have no recognizable one.
- An ordination is a view. Subfamily assignment stays with the workflow, and any group seen in a plot has to be shown as an alignment before it means anything.

### 4.1 Repair the toy

Cheap, and needed before the table in 2.3 can be quoted. It needs an environment with `umap-learn` and `pacmap`.

| Change | Why | Known answer to check against |
|---|---|---|
| Score in sequence space: cluster each picture, then compute silhouette of those clusters on the K2P matrix | removes the circularity described in 2.5 | silhouette of the true labels on K2P is one fixed number, the same for every method; no method can honestly exceed what that allows |
| Give every method the same input (all K2P, then all PCA-30) | separates method from input | PCA on one-hot must reproduce 0.905 |
| Run NMDS from the PCoA solution, with several starts and more iterations | the present NMDS result may be a failed run | stress must not rise above 0.251 |
| Repeat over several seeds and report the spread | one seed cannot rank methods that differ by 0.01 | seed 1 must reproduce the present table |
| Print the negative-eigenvalue share with every PCoA percentage | 3.3 | for the toy Jaccard matrix it must be 0 |
| Treat missing cells as missing in distances and counts | at present they are counted as absent | with 0 % missing the results must not change |
| Compare the UPGMA tree with the simulated tree by a tree distance | "recovers clades" is by eye | the true tree is in the script |

### 4.2 Make the toy as hard as real copies

Add indels, 5′-truncated copies, faster decay at CpG sites, and subfamilies of very unequal size. Then repeat 4.1.

- Question: do the non-linear methods still separate subfamilies when copies are gapped and truncated, and does K2P (which ignores gapped columns) remain a usable distance?
- Check: ARI has to fall below the present 0.99. If it does not, the toy is still too easy to rank anything.

### 4.3 One real family from a SINEderella run

- Input: the copies of one family that the workflow accepted (unique flanks), in the alignment the workflow already makes. Nothing new is searched or aligned.
- Picture: PCA and one non-linear method side by side, colored by the subfamily the workflow assigned.
- Question: do the assigned subfamilies appear as groups, and which copies sit between groups?
- Control: a family whose subfamilies are already established should show them. If it does not, the picture is not informative for this kind of data and the idea stops here.
- Open point to decide first: how truncated copies enter the distance (pairwise over shared columns, or full-length copies only).

### 4.4 A real presence/absence matrix

The [sq2 pilot](../../projects/sq2-distant-orthologous-loci/) is a plan; it has produced no matrix. The tree + heatmap and the UpSet (old next step 2) wait for it.

A shallower matrix may be closer: the [TE-insertion anchor graph](../../projects/te-anchor-graph-pansineome/) is implemented and has run on *Darevskia* (149,234 groups). If its group table records presence per genome, it is already a presence/absence matrix, and the heatmap script could be pointed at it. This has not been checked.

In either case the matrix keeps three states (present, empty site, missing), and missing is never turned into absent.

### 4.5 Insertion phylogeny from the matrix

Only after 4.4. The methods exist and have been used on SINE data:

- **Dollo parsimony** (PHYLIP `dollop`), as used for 10,000+ Ves insertions in *Myotis* (Korstian et al. 2022), together with SplitsTree networks.
- **Binary maximum likelihood** in IQ-TREE: models JC2 and GTR2, with the +ASC correction because a locus table contains no constant sites. IQ-TREE has no Dollo model.
- **KKSC test** (Kuritzin et al. 2016) for nodes where insertions conflict: it separates incomplete lineage sorting from ancestral hybridization using the counts of conflicting insertions.

Control: nodes that are not in dispute must be recovered. If they are not, the absence calls are wrong and the tree means nothing.

### 4.6 Locus strip in the viewer

Port the raster locus strip (figure 6, left panel) to the MSA / SINE viewer with a canvas (old next step 3). The real version draws real flank matches, so it also shows at a glance whether an "empty site" row has continuous flanks.

### 4.7 Looked at and set aside

These were raised in an outside review of this page.

| Idea | Why not now |
|---|---|
| Sequence graphs (pggb) across distant species | the pggb README gives examples up to about 10 % divergence; flanks between distant squamates are far beyond that, and it would need long alignments. The anchor graph already covers the shallow case with short flanks |
| GraffiTE | built for insertion polymorphism inside a species (shown on human, *Drosophila melanogaster*, maize, *Cannabis sativa*) |
| Language-model embeddings of copies | nothing to gain on a toy that is already solved; the published deep-learning TE classifiers work at order / superfamily level, not at subfamily level. Worth a look only if 4.3 fails |
| Requiring a TSD for every "present" call | would drop true old insertions; see the rules above |
| SVDquartets on the presence/absence matrix | it is a method for SNP data |
| Model-based inference from Heaps' α | α is descriptive (3.3) |
| Somatic insertions in cancer; regulatory exaptation | different questions, different data; not part of this page |

## 5. Literature

Checked on 2026-11-11 by opening each source. "Abstract" means only the abstract was read.

**Ordination**

- Kobak D, Berens P (2019) The art of using t-SNE for single-cell transcriptomics. *Nat Commun* 10:5416. [doi:10.1038/s41467-019-13056-x](https://doi.org/10.1038/s41467-019-13056-x). Full text: "the relative position of clusters on the t-SNE plot is almost arbitrary", depending mostly on the random seed.
- Kobak D, Linderman GC (2021) Initialization is critical for preserving global data structure in both t-SNE and UMAP. *Nat Biotechnol* 39:156-157. [doi:10.1038/s41587-020-00809-z](https://doi.org/10.1038/s41587-020-00809-z). Preprint text read.
- Chari T, Pachter L (2023) The specious art of single-cell genomics. *PLoS Comput Biol* 19:e1011288. [doi:10.1371/journal.pcbi.1011288](https://doi.org/10.1371/journal.pcbi.1011288). Reduction to 2 dimensions "inevitably induces significant distortion".
- Wang Y, Huang H, Rudin C, Shaposhnik Y (2021) Understanding how dimension reduction tools work. *JMLR* 22(201):1-73. [jmlr.org/papers/v22/20-1061.html](https://jmlr.org/papers/v22/20-1061.html). Abstract. Source of PaCMAP.
- Negative eigenvalues in PCoA and the Lingoes / Cailliez corrections: documentation of [`ape::pcoa`](https://search.r-project.org/CRAN/refmans/ape/html/pcoa.html).
- NMDS stress above 0.2 as "poor": attributed to Kruskal (1964, *Psychometrika*). Seen only second-hand, in a software manual; the original was not opened.

**Subfamily clustering**

- Orozco-Arias S, Sierra P, Durbin R, González J (2024) MCHelper automatically curates transposable element libraries across eukaryotic species. *Genome Res* 34:2256-2268. [doi:10.1101/gr.278821.123](https://doi.org/10.1101/gr.278821.123). Splits a family into subfamilies by clustering on K2P distances from the alignment; the code uses DBSCAN.
- Shingleton AM (2022) Subfamily Clustering Using Label Uncertainty (SCULU). MS thesis, University of Montana. [scholarworks.umt.edu/etd/11913](https://scholarworks.umt.edu/etd/11913). It merges subfamilies that cannot be told apart in annotation; it does not find new ones.
- TEclass (Abrusán et al. 2009, *Bioinformatics* 25:1329), DeepTE (Yan et al. 2020, *Bioinformatics* 36:4269), NeuralTE and BERTE (2024, preprints), TEclass2 (Bickmann et al. 2025, *Biology* 15:59): all classify at class, order or superfamily level. Abstracts.

**Insertion presence/absence as phylogenetic markers**

- Korstian JM, Paulat NS, Platt RN, Stevens RD, Ray DA (2022) SINE-based phylogenomics reveal extensive introgression and incomplete lineage sorting in *Myotis*. *Genes* 13:399. [doi:10.3390/genes13030399](https://doi.org/10.3390/genes13030399). Full text.
- Kuritzin A, Kischka T, Schmitz J, Churakov G (2016) Incomplete lineage sorting and hybridization statistics for large-scale retroposon insertion data. *PLoS Comput Biol* 12:e1004812. [doi:10.1371/journal.pcbi.1004812](https://doi.org/10.1371/journal.pcbi.1004812). Abstract.
- Grechko VV, Kosushkin SA, Borodulina OR, Butaeva FG, Darevsky IS (2011) Short interspersed elements (SINEs) of squamate reptiles (Squam1 and Squam2): structure and phylogenetic significance. *J Exp Zool B* 316B:212-226. [doi:10.1002/jez.b.21391](https://doi.org/10.1002/jez.b.21391).
- IQ-TREE binary models and +ASC: [iqtree.github.io/doc/Substitution-Models](https://iqtree.github.io/doc/Substitution-Models).
- Chifman J, Kubatko L (2014) Quartet inference from SNP data under the coalescent model. *Bioinformatics* 30:3317-3324. [doi:10.1093/bioinformatics/btu530](https://doi.org/10.1093/bioinformatics/btu530). Abstract.

**Pangenome visuals and graphs**

- Panache: Durant É, Sabot F, Conte M, Rouard M (2021) *Bioinformatics* 37:4556-4558. [doi:10.1093/bioinformatics/btab688](https://doi.org/10.1093/bioinformatics/btab688). Browser-style viewer of presence/absence along a linearized pangenome. Abstract.
- PanViz: Pedersen TL, Nookaew I, Ussery DW, Månsson M (2017) *Bioinformatics* 33:1081-1082. [doi:10.1093/bioinformatics/btw761](https://doi.org/10.1093/bioinformatics/btw761). Abstract.
- Heaps' law for open / closed pangenomes: attributed to Tettelin H, Riley D, Cattuto C, Medini D (2008) *Curr Opin Microbiol* 11:472-477. Not opened. Tettelin et al. 2005 (*PNAS* 102:13950), often cited for it, fits an exponential and does not use Heaps' law.
- pggb: Garrison E, Guarracino A, et al. (2024) Building pangenome graphs. *Nat Methods* 21:2008-2012. [doi:10.1038/s41592-024-02430-3](https://doi.org/10.1038/s41592-024-02430-3). Divergence examples are from the [README](https://github.com/pangenome/pggb/blob/master/README.md).
- GraffiTE: Groza C, Chen X, Wheeler TJ, Bourque G, Goubert C (2024) *Nat Commun* 15:8915. [doi:10.1038/s41467-024-53294-2](https://doi.org/10.1038/s41467-024-53294-2). Abstract.

## 6. Reproduce

Run the two scripts in `projects/sine-ordination-toy/`; each writes its figures to the working directory. Versions used: Python 3.13, numpy, scipy, scikit-learn, matplotlib, umap-learn 0.5.12, pacmap 0.9.1. `ordination_toy.py` also writes `ordination_metrics.json`; `pansineome_toy.py` writes `toy_pansineome_PAV.tsv`.
