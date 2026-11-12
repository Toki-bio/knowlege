---
title: Beyond PCA for SINEderella output; pan-SINEome visuals
type: notes
---

> Type: notes

# Beyond PCA for SINEderella output, and visuals for orthologous loci / pan-SINEome

**Status:** two toy ordination simulations and one real family (*Scalopus aquaticus* Tal SINE, 2.8) for SINE copies; toy data only for presence/absence of loci.

Code, simulated matrix and metrics: [`projects/sine-ordination-toy/`](https://github.com/Toki-bio/knowlege/tree/main/projects/sine-ordination-toy).

Related: [dimensionality reduction for SINE separation](../dim-reduction-for-sine-separation/), [clustering for SINE subfamily discovery](../clustering-for-sine-subfamily-discovery/), [HDBSCAN vs vsearch for SINE grouping](../hdbscan-vs-vsearch-for-sine-grouping/), [TE-insertion anchor graph](../../projects/te-anchor-graph-pansineome/), [sq2 distant orthologous loci](../../projects/sq2-distant-orthologous-loci/).

## Conclusion: how to present dimensionality reduction of SINE copies

Based on the toys (2.6, 2.7, 2.9) and four real families (2.8–2.10). An earlier version of this section used K2P over all columns as the only input; 2.9 shows that was too narrow.

1. **The input matters more than the picture method.** A subfamily is a few diagnostic positions among many private mutations. A distance over all columns (K2P) dilutes them. Two inputs should be made and shown side by side:
   - **all columns:** K2P with pairwise deletion, columns covered by fewer than half the copies dropped, CpG columns masked;
   - **diagnostic positions:** minor-allele characters (substitutions, and indels as characters) that co-vary with a character at another position, found without labels (2.9), then a Hamming distance on them with pairwise deletion.

   Not one-hot, and not PCA of everything: truncated copies are then placed by their missing block, and a signal carried by a few columns stays below the noise of the rest.
2. **Picture:** UMAP (`n_neighbors` 30, `min_dist` 0.3) with t-SNE (perplexity 40) beside it, of each input. Colours = the workflow's subfamily assignment. Axes without units or percentages.
3. **Numbers under the picture, measured on the distance and not on the 2D picture:**
   - how many positions were selected as co-varying (none = no evidence of subfamilies in this set);
   - the share of copies whose 15 nearest copies are mostly of their own assigned subfamily;
   - the silhouette of the assigned subfamilies.
4. **How to read it:**
   - The diagnostic-position distance is the better measure of whether subfamilies are distinct: in the toy it equals a distance on the true diagnostic columns. Its picture is not better than the K2P picture, because a few discrete positions break a group into islands; do not count islands.
   - If the all-column picture shows two blobs and the diagnostic-position numbers show more, as for saq, the subfamilies are real at the level of positions but differ at too few of them for an all-column distance.
   - Groups proposed from a picture or from the nested procedure of 2.9 are candidates. In the toy the nested procedure returned one or two groups too many. Each candidate has to be shown as an alignment before it means anything.
   - A picture never assigns copies and never defines a subfamily.

Limits: not tested with recombination between subfamilies or very unequal subfamily sizes; four real families (2.10), 100 copies per subfamily; the nested procedure over-splits in the toy for a reason not yet found.

## In short

- A 2D PCA of SINE copies shows a low percentage on its axes because the subfamily signal is spread over many weak axes. The low percentage is expected, and by itself it does not make the PCA picture wrong.
- On the toy, t-SNE, UMAP and PaCMAP recovered the six simulated subfamilies in each of 5 seeds (ARI 0.98 to 0.99, sd at most 0.006). PCA, PCoA and NMDS did so far less reliably (ARI 0.73 ± 0.18, 0.69 ± 0.15 and 0.22 ± 0.10). The first single-seed table in 2.3 (PCA ARI 0.905) happened to use a favourable seed. See [2.6](#26-repaired-run-5-seeds).
- Scored in sequence space, t-SNE, UMAP and PaCMAP all sit at the silhouette value of the true labels themselves (about 0.10). They reproduce the true grouping and cannot exceed it, and the three cannot be ranked against each other. The much higher silhouettes in the first table (0.86 to 0.89) were measured in the 2D picture.
- With indels, 5′ truncation and CpG decay added (2.7), only t-SNE and UMAP on K2P distances with pairwise deletion stayed good (ARI 0.93). Anything built on one-hot fell to 0.30–0.46, because truncated copies are placed by their missing block. Masking CpG columns improved every K2P picture. Without being told the number of subfamilies, HDBSCAN on the UMAP picture found the 6 subfamilies in 13 of 15 runs; HDBSCAN on the distances themselves found almost nothing. See [2.7](#27-harder-toy-indels-truncation-cpg-decay-number-of-subfamilies-not-given).
- A distance on positions that co-vary across copies (found without labels) is far better than K2P over all columns: on saq the K2P picture shows 2 groups, while the co-varying positions recover most of the published s1–s8. See [2.9](#29-other-inputs-diagnostic-positions-found-by-covariation).
- Not yet in any toy: recombination between subfamilies, very unequal subfamily sizes.
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

### 2.3 Results of the first run (seed 1 only)

Superseded by [2.6](#26-repaired-run-5-seeds); kept because the figures below come from this run.

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
- **PCA separates most of it in a good seed** (ARI 0.905 in seed 1, but 0.733 ± 0.178 over 5 seeds), with the neighboring subfamilies overlapping.
- **t-SNE, UMAP and PaCMAP separate all six subfamilies** (ARI 0.981 to 0.993 in seed 1; 0.98 to 0.99 in every seed, see 2.6).
- **Trustworthiness is 0.67 to 0.71 for every method**, so it did not tell the methods apart here.
- UMAP, t-SNE and PaCMAP have no variance percentage, and their axes have no meaning. The position of one group relative to another should not be read.

### 2.5 What the toy cannot show

- **Three of the four scores in the first table are measured in the 2D picture (repaired in 2.6).** Silhouette, kNN accuracy and ARI all reward a method for drawing compact, separate groups. That is what UMAP and PaCMAP are built to do, so their high silhouette (0.86, 0.89) is partly a property of the method. Only trustworthiness refers back to sequence distances, and it does not separate the methods.
- **The methods did not get the same input.** PCA used one-hot, PaCMAP used the first 30 PCs, the others used K2P. Method and input cannot be told apart in this table.
- **NMDS was run weakly.** One random start and 300 iterations; the NMDS panel is star-shaped, which looks like a run that did not converge. "NMDS did not beat PCA" may describe this run and not the method. Checked in 2.6: restarted from the PCoA solution with 3,000 iterations, NMDS did not improve (ARI 0.22 ± 0.10), so the weak run does not explain it.
- **Initialization differs.** t-SNE started from random positions and UMAP from its default. The large-scale arrangement of groups depends on initialization (Kobak & Linderman 2021), so the arrangements in the figure are not comparable.
- **The toy is easy.** No indels, no truncated copies, one noise model, one seed. An ARI of 0.99 leaves no room to compare anything better.
- **K is given.** ARI uses Ward clustering with the true number of subfamilies (6). On real data that number is unknown.

### 2.6 Repaired run, 5 seeds

Script `ordination_toy_v2.py`, run on therioserver on 2026-11-11 (Python 3.14, numpy 2.5.4, scikit-learn 1.9.1, umap-learn 0.5.12, pacmap 0.9.1). Results: `ordination_v2_metrics.json`; log: `ordination_v2_run_seed5.log`.

What changed against the first run:

- **Scored in sequence space.** Each picture is clustered (Ward, k = 6) and the silhouette of those clusters is computed on the K2P matrix. For reference, the silhouette of the *true* labels on K2P is 0.098 (mean over seeds).
- **Same input for every method.** Condition O: PCA-30 of the one-hot matrix (PCA itself on one-hot). Condition K: the K2P matrix (PaCMAP gets the first 30 PCoA axes of K2P).
- **NMDS** starts from the PCoA solution (`init="classical_mds"`), 3,000 iterations.
- **t-SNE** also run from the PCoA start.
- **Seeds 1 to 5.** The seed sets both the simulation (including which sites are diagnostic) and the method; the table gives mean ± sd (population sd, n = 5).
- **Control passed.** Seed 1 reproduces PCA ARI 0.905 and PC1 + PC2 = 11.9 %.

| Picture | Trustworthiness | Silhouette in K2P | kNN accuracy | ARI |
|---|---|---|---|---|
| true labels (reference) | | 0.098 | | |
| O: PCA (one-hot) | 0.672 ± 0.018 | 0.030 ± 0.028 | 0.930 ± 0.052 | 0.733 ± 0.178 |
| O: t-SNE | 0.711 ± 0.004 | 0.101 ± 0.007 | 0.992 ± 0.001 | 0.980 ± 0.002 |
| O: UMAP | 0.680 ± 0.004 | 0.101 ± 0.008 | 0.992 ± 0.002 | 0.981 ± 0.004 |
| O: PaCMAP | 0.673 ± 0.004 | 0.101 ± 0.005 | 0.993 ± 0.002 | 0.982 ± 0.004 |
| K: PCoA | 0.683 ± 0.022 | 0.053 ± 0.033 | 0.877 ± 0.078 | 0.692 ± 0.149 |
| K: NMDS (PCoA start) | 0.678 ± 0.008 | 0.080 ± 0.026 | 0.912 ± 0.032 | 0.224 ± 0.096 |
| K: t-SNE | 0.696 ± 0.004 | 0.100 ± 0.003 | 0.996 ± 0.001 | 0.989 ± 0.004 |
| K: t-SNE (PCoA start) | 0.697 ± 0.004 | 0.099 ± 0.004 | 0.996 ± 0.001 | 0.990 ± 0.002 |
| K: UMAP | 0.673 ± 0.003 | 0.100 ± 0.004 | 0.996 ± 0.001 | 0.991 ± 0.002 |
| K: PaCMAP | 0.691 ± 0.004 | 0.096 ± 0.006 | 0.992 ± 0.002 | 0.980 ± 0.006 |

Reading:

- **t-SNE, UMAP and PaCMAP are equivalent here.** ARI 0.98 to 0.99 in every seed, and the differences between them are within the sd. Their silhouette in K2P (0.096 to 0.101) equals that of the true labels: they recover the true grouping and cannot exceed it.
- **The silhouette in the picture overstated them.** In the first table UMAP and PaCMAP had 0.86 and 0.89; that measures compact blobs in 2D, not agreement with sequence distances.
- **PCA, PCoA and NMDS depend strongly on the seed.** PCA ARI per seed: 0.905, 0.920, 0.441, 0.761, 0.641. PCoA: 0.592, 0.821, 0.447, 0.777, 0.822. NMDS: 0.136, 0.082, 0.303, 0.282, 0.316. The original single-seed PCA value (0.905) was close to the best of five.
- **Input choice hardly matters.** The same method on K2P or on PCA-30 of one-hot differs by about 0.01 in ARI. Starting t-SNE from the PCoA solution instead of random changes ARI by 0.001.
- **Silhouette in K2P is low for every method** (at most 0.11), because the 9 % per-copy noise is much larger than the 6 diagnostic sites. It cannot separate the good methods from each other; ARI separates them from PCA.
- **Seed 3 is a harder draw.** PCA and PCoA fall to 0.44 ARI, the true-label silhouette is lowest (0.087) and the negative-eigenvalue share highest (0.378). That the subfamilies lie closer together in that draw is a guess; the cause was not checked.
- **NMDS stress in seed 3 is 0.022, against 0.245 to 0.251 in the other four seeds,** while its ARI stays poor (0.303). That run is not understood and its stress should not be trusted; a degenerate solution is a guess, not a finding.
- **Negative-eigenvalue share of the K2P PCoA** is 0.331 to 0.378 in all seeds.

Still open: no indels, no truncated copies, one noise model; the number of subfamilies is given to the clustering; ARI measures whether Ward clustering of the picture returns the true groups, which is a property of the picture as much as of the method.

### 2.7 Harder toy: indels, truncation, CpG decay; number of subfamilies not given

Script `ordination_toy_v3.py`, therioserver, 2026-11-11, same versions as 2.6. Results: `ordination_v3_base.json`, `ordination_v3_notrunc.json`, `ordination_v3_maskcpg.json` and the matching `.log` files.

**Why this toy.** I found no published SINE subfamily study that tested its method on simulated copies with known labels; they validate against an established scheme, trees, copy age, PCR of young loci or function (literature in 5). The only per-copy ordination case found, mouse L1 promoter monomers (Zhou et al. 2019), excluded truncated copies, used PCA to 100 axes, then HDBSCAN without a set number of clusters, and t-SNE only for display. So the toy borrows its ingredients from TE simulators and mutation studies:

| Ingredient | Value | Source |
|---|---|---|
| Subfamily tree, sizes, per-copy divergence | as in 2.1 | |
| CpG decay | substitution rate ×6 at CpG sites of the copy's consensus, as C→T or G→A | Xing et al. 2004 (CpG : non-CpG ≈ 6 in AluY/AluS) |
| Transitions | 2/3 of other substitutions | my choice |
| Indels | 0.15 events per substitution; 40 % insertions, 60 % deletions; length geometric, mean ≈ 1.7, max 20 | between Lunter et al. 2006 (≈ 1:10) and Hubley et al. 2022 (0.20; 0.08 : 0.12; mean 1.7, max 20) |
| 5′ truncation | 30 % of copies, 0–50 % of the length removed | TEvarSim defaults (Miao & Li 2026) |

Deletions and truncations become gap columns. Insertions add no column, so they never enter a distance (a simplification). Measured in the base version: 8 % of cells are gaps, 27–32 % of copies are truncated, CpG sites changed 6.0 times as often as other sites (local check, seeds 1 and 3). Distance: K2P with pairwise deletion (only columns where both copies have a base; fewer than 50 shared columns counts as saturated). One-hot: a gap is all zeros.

Three versions, 5 seeds each: **base**; **no truncation**; **CpG masked** (columns at CpG sites of the family consensus removed, as COSEG users are told to discount TG/CA changes). Control: the base version reproduced, to the third decimal, a separate earlier run of the same script.

ARI (mean ± sd over 5 seeds; Ward clustering, number of subfamilies given):

| Picture | base | no truncation | CpG masked |
|---|---|---|---|
| one-hot: PCA | 0.30 ± 0.05 | 0.30 ± 0.06 | 0.37 ± 0.01 |
| one-hot PCA-30: UMAP | 0.46 ± 0.05 | 0.91 ± 0.06 | 0.52 ± 0.06 |
| one-hot PCA-30: PaCMAP | 0.44 ± 0.07 | 0.90 ± 0.08 | 0.48 ± 0.03 |
| K2P: PCoA | 0.32 ± 0.06 | 0.38 ± 0.07 | 0.65 ± 0.18 |
| K2P: NMDS | 0.26 ± 0.04 | 0.28 ± 0.02 | 0.32 ± 0.05 |
| K2P: t-SNE | **0.93 ± 0.04** | 0.98 ± 0.02 | **0.96 ± 0.04** |
| K2P: UMAP | **0.93 ± 0.05** | 0.99 ± 0.01 | **0.95 ± 0.05** |
| K2P: PaCMAP (PCoA-30) | 0.78 ± 0.12 | 0.80 ± 0.11 | **0.95 ± 0.04** |
| true-label silhouette on K2P | 0.115 | 0.111 | 0.155 |

HDBSCAN, number of subfamilies **not** given (`min_cluster_size` 20; unclustered copies counted as one extra label):

| Clustered | base | no truncation | CpG masked |
|---|---|---|---|
| UMAP picture of K2P | 0.90 ± 0.09; 6 clusters in 4 seeds, 5 in 1; 3 % unclustered | 0.99 ± 0.01; 6 in all; 0 % | 0.93 ± 0.09; 6, 6, 5, 6, 6; 1 % |
| t-SNE picture of K2P | 0.83 ± 0.10; 6, 6, 5, 6, 6; 8 % | 0.96 ± 0.04; 6 in all; 2 % | 0.91 ± 0.09; 6, 6, 5, 6, 6; 3 % |
| K2P matrix itself | 0.00; 0–2 clusters; 99 % unclustered | 0.00; none; 100 % | 0.00; 2–3; 88 % |
| one-hot PCA-100 (Zhou et al. 2019 recipe) | 0.00; 2 clusters; 67 % | 0.00; 2–5; 91 % | 0.01; 2–6; 57 % |

What it shows:

- **One-hot fails with truncated copies.** UMAP and PaCMAP on one-hot drop from 0.90 (no truncation) to 0.45: a truncated copy has a block of zeros and is placed by that block, not by its subfamily. Zhou et al. removed truncated copies before PCA for this reason. K2P with pairwise deletion compares only shared columns and is affected far less (0.93 against 0.99).
- **On K2P, t-SNE and UMAP stay best** (0.93 in base). PaCMAP is weaker and less stable (0.78 ± 0.12), unlike in the easy toy.
- **Masking CpG columns helps every K2P picture** (PCoA 0.32 → 0.65, PaCMAP 0.78 → 0.95) and raises the true-label silhouette from 0.115 to 0.155: CpG decay is noise shared across subfamilies.
- **PCA, PCoA and NMDS do not separate the subfamilies in this toy** (ARI 0.26 to 0.32 in base).
- **Without the number of subfamilies, HDBSCAN on the UMAP picture found 6 groups in 13 of 15 runs and 5 in 2** (all in seed 3), with on average 0 to 3 % of copies left unclustered (up to 10 % in a single run). HDBSCAN directly on the K2P matrix or on 100 PCA axes found almost nothing, probably because in that many dimensions the density contrast between subfamilies, which differ at 6 sites among 9 % noise, is too weak. Clustering the 2D picture works here, but the picture's clusters are partly made by UMAP itself (see 2.6), so on real data every such group has to be shown as an alignment before it means anything.
- **Seed 3 is the hardest draw in the base and CpG-masked versions** (lowest UMAP ARI in both, lowest t-SNE ARI in base; 5 instead of 6 clusters in both), but not without truncation, where it is among the best. Not explained.

Seed 3 of the easy toy (2.6) has a separate, explained cause: one pair of copies whose K2P distance is saturated (the formula's logarithm hits its clip value). Removing them alone moves seed 3 back in line: PCA ARI 0.441 → 0.917, PCoA 0.447 → 0.794, NMDS stress 0.022 → 0.251, true-label silhouette 0.087 → 0.092. So two extreme copies were enough to distort PCA, PCoA and NMDS; t-SNE, UMAP and PaCMAP were unaffected. The low NMDS stress was a by-product of those outliers, not a good fit. (In the harder toy no copy is saturated, so its "filter" arm equals the base and is not shown.)

Still not in any toy: gene conversion and recombination between subfamilies, copies that continue to decay in a continuum rather than forming groups, and very unequal subfamily sizes. Carey et al. 2021 discuss recombination and chance mutation at diagnostic sites as reasons subfamily labels disagree on real data; Wacholder et al. 2014 found recombination in 7.7 % of full-length LAVA copies.

### 2.8 First real family: *Scalopus aquaticus* Tal SINE, published subfamilies s1–s9

Script and results: `projects/sine-ordination-toy/real_saq/` (`saq_ordination.py`, `saq_ordination.json`, `saq_ordination.png`, `saq_hdbscan_vs_labels.tsv`, `saq_between_copies.tsv`, `saq_all.aln.fa`). Therioserver, 2026-11-11.

**Input.** The 100 random copies per subfamily published in the Tal repository (`saq/alignments/saq_s*_rand100.aln.fa`, last changed 2026-04-28). These rows carry 50 bp of left and 70 bp of right flank; the element body was cut out by that rule for the 866 copies whose length is exactly element + 120 bp (34 copies at contig ends dropped). The flanks of these copies were **not** tested for uniqueness in this run. Bodies + 9 consensuses aligned with MAFFT 7.526 `--auto`. Distance as in 2.7: K2P with pairwise deletion, columns with a base in fewer than half the copies dropped (251 kept), CpG columns of the majority consensus masked (30), 221 columns used.

**What came out.**

- The published labels do not form groups in sequence distance: their silhouette on K2P is 0.01 (toy reference values 0.10–0.16). Ward clustering into 9 on the pictures: ARI 0.28 (PCoA), 0.45 (t-SNE), 0.37 (UMAP).
- HDBSCAN on the UMAP picture, number not given, found **2 clusters** with no copy left out. Cluster 0: all of s1, s2, s7g, s9 and 64 of 99 s3. Cluster 1: all of s4, s5, s6, s8 and 35 of s3.
- The 9 consensuses fall into the same two sets. Within {s1, s2, s3, s7g, s9}: 4–18 substitutions apart; within {s4, s5, s6, s8}: 5–16; between the sets: 20–34 substitutions and 29–53 columns differing by an indel.
- Mean K2P within one subfamily is 0.07–0.23, far more than the distance between neighbouring consensuses (for example s1 and s9 differ at 4 sites). So at the copy level the close subfamilies overlap, as in the toy when noise swamps diagnostic sites.
- Nearest consensus (mismatches plus gap differences) agrees with the published label for 677 of 866 copies. Largest disagreements: s3 43 of 99 copies nearest s8 (the other set); s9 45 of 99 nearest its own, the rest nearest s1, s7g or s2; s6 24 of 94 nearest s4; s5 18 of 90 nearest s4.
- In the t-SNE and UMAP pictures the s3 copies sit in two places, one beside each set.

**Labels.** The s1–s9 labels of these rows come from the run's automatic assignment of copies to consensuses by bit score, not from hand curation.

**Flank check (within the 866 copies).** Method of the Tal `eri/LOG.md` test: position-by-position flank identity between every pair, no alignment. No pair shares flanks above 60 % on both sides, also when shifts of up to ±5 bp are allowed (`saq_flank_test.py`). A copy whose flanks match a locus outside this sample would not be seen; that genome-wide test needs the saq genome (on KIT) and was not run.

**s3 alignment:** the 99 s3 copies in their published order, with the s3, s8 and s7g consensuses as the first rows; each name carries its nearest consensus and its distances to s3 and s8 (`saq_s3_copies_vs_s3_s8_s7g.aln.fa`, columns from the joint MAFFT alignment). [Open in MSA viewer](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2Fknowlege%2Fmain%2Fprojects%2Fsine-ordination-toy%2Freal_saq%2Fsaq_s3_copies_vs_s3_s8_s7g.aln.fa&title=saq%20s3%20copies%20vs%20s3%2C%20s8%2C%20s7g%20consensuses).

**Owner's reading of the s3 alignment (2026-11-12):** s3 and s8 split faintly; a signal from s7g is not obvious; the separation is very weak.

The low scores in this section belong to the input, K2P over all columns, and not to the family: see 2.9.

### 2.9 Other inputs: diagnostic positions found by covariation

Owner's point (2026-11-12): the input should not be limited to a K2P distance; fine-scale structure and a matrix with heavy diagnostic positions should be considered. Script `sine_inputs_v4.py`; saq: `real_saq/saq_inputs_v4.py`. Therioserver, 2026-11-12. Results: `inputs_v4_d6c`, `_d3c`, `_d2c`, `_nosubc` (`.json`, `.log`); `real_saq/saq_inputs_v4.json`, `saq_inputs_v4_features.tsv`, `saq_inputs_v4_crosstabs.txt`, `saq_inputs_v4.png`.

**Idea.** The diagnostic positions of a subfamily co-vary across copies: a copy that has one usually has the others. Private mutations do not. So positions can be selected without labels.

**Method.**

- Characters: every (column, minor state) present in at least 3 % of the copies. Terminal gaps are "not observed"; internal gaps and insertions can be used as characters, a run of adjacent gap columns counting once.
- Covariation: correlation between two characters at different positions (at least 3 columns apart, 25 if an indel is involved), over copies where both are observed.
- Copy age: an old copy carries more of everything, which alone makes all noise positions co-vary. So each character is first centred within bins of copies of similar divergence from the majority consensus. The threshold is the largest value seen after shuffling each character within those bins.
- A character is **selected** if it co-varies above the threshold with at least one other.
- **Nested:** copies are split on the selected characters (HDBSCAN on UMAP proposes parts; parts are merged until every two differ at 2 or more positions with a frequency difference of at least 0.5); then selection is repeated inside each part, until nothing co-varies.

**Toy** (the harder toy of 2.7 with 6, 3 or 2 diagnostic sites per step of the tree, and a control with no subfamilies; 5 seeds each):

| | 6 sites | 3 sites | 2 sites | no subfamilies |
|---|---|---|---|---|
| true diagnostic columns found, per seed | 30/30, 27/27, 28/28, 29/29, 28/28 | 15/15 ×4, 14/14 | 10/10 ×4, 9/9 | none exist |
| other columns selected, per seed | 7, 8, 8, 7, 9 | 10, 1, 4, 4, 5 | 3, 1, 3, 2, 2 | 0, 0, 2, 0, 0 |

In the one seed examined (3 sites, seed 1), the other columns were mostly the CpG partner of a diagnostic site: the diagnostic change made or destroyed a CpG, so the neighbour mutates at a different rate in that subfamily.

Scores on the distance itself (mean of 5 seeds): share of copies whose 15 nearest copies are mostly their own subfamily / silhouette of the true labels.

| Input | 6 sites | 3 sites | 2 sites |
|---|---|---|---|
| K2P | 0.98 / 0.12 | 0.90 / 0.05 | 0.84 / 0.04 |
| K2P, CpG masked | 0.98 / 0.16 | 0.91 / 0.08 | 0.87 / 0.06 |
| all characters, unweighted | 0.97 / 0.13 | 0.90 / 0.07 | 0.84 / 0.05 |
| selected characters | 0.98 / 0.44 | 0.94 / 0.51 | 0.91 / 0.51 |
| all characters, selected ones carrying 75 % of the weight | 0.98 / 0.38 | 0.95 / 0.42 | 0.91 / 0.42 |
| true diagnostic columns (upper bound) | 0.99 / 0.50 | 0.95 / 0.52 | 0.94 / 0.51 |
| PCA of all characters, axes by parallel analysis | 0.77 / 0.07 | 0.20 / −0.09 | 0.19 / −0.08 |

Groups recovered (ARI against the truth; 6 true subfamilies):

| Route | 6 sites | 3 sites | 2 sites |
|---|---|---|---|
| UMAP of K2P CpG-masked, Ward k = 6 | 0.95 | 0.77 | 0.69 |
| UMAP of selected characters, Ward k = 6 | 0.67 | 0.68 | 0.58 |
| UMAP of all characters, selected ones carrying 50 % of the weight, Ward k = 6 | 0.93 | 0.72 | 0.58 |
| UMAP of the true diagnostic columns, Ward k = 6 | 0.93 | 0.75 | 0.69 |
| HDBSCAN on UMAP of K2P CpG-masked, number not given | 0.93 (6, 6, 5, 6, 6 groups) | 0.60 (5, 5, 3, 2, 5) | 0.52 (5, 3, 4, 4, 5) |
| nested, number not given | 0.91 (8, 7, 8, 8, 7) | 0.86 (6, 6, 8, 8, 7) | 0.72 (7, 8, 6, 5, 6) |

Nested on the no-subfamily control: one group in 5 of 5 seeds.

What the toy shows:

- **Selection works.** Nearly every true diagnostic column is found, with few others, down to 2 sites per step. With no subfamilies almost nothing is selected and the nested procedure returns one group.
- **The selected distance is as good as knowing the diagnostic columns** (silhouette 0.44–0.51 against 0.50–0.52; K2P 0.04–0.16).
- **The picture does not follow.** A UMAP of 10–40 discrete positions breaks subfamilies into islands, and even the true-column picture is no better than the K2P picture. Keeping all characters with the selected ones heavy repairs the strong-signal case (0.93) and not the weak ones. The gain of the selected input is in the distance and in the nested grouping, not in the 2D picture.
- **Nested grouping is much better than clustering the K2P picture when the signal is weak** (0.86 against 0.60; 0.72 against 0.52), **and returns one or two groups too many when it is strong.** A stricter merging rule (the separating positions must co-vary with each other in the parent set) changed one run of fifteen. The cause is not found.
- **PCA of all characters fails**: a signal in a few columns does not rise above the noise of the others.
- Weighting the selected characters by the strength of their covariation was worse than weighting them equally.

**saq, 866 copies** (same alignment as 2.8; label-free unless said):

- Selected: 79 of 207 substitution characters; with indels 168 of 307 (90 of them indels).
- **Check against the published consensuses, which the selection never saw:** of the 87 columns where the 9 consensuses differ and a character exists, 80 were selected. In the other direction, 78 of the 168 selected characters are a state in which the consensuses differ; 54 selected columns are ones where they do not differ. What those 54 are (further structure, CpG partners, alignment shifts next to indels) was not examined.
- **Split-half:** selection repeated in two disjoint random halves of the copies, 5 times: 78–99 characters selected in both halves, of 89–139 per half.

| Input | 15 nearest mostly own published subfamily | silhouette of published labels | HDBSCAN groups on UMAP |
|---|---|---|---|
| K2P, CpG masked (2.8) | 0.64 | 0.01 | 2 |
| all characters | 0.72 | 0.03 | 2 |
| selected, substitutions | 0.81 | 0.09 | 2 |
| selected, substitutions + indels | 0.83 | 0.07 | 3 |
| columns where the published consensuses differ (not label-free) | 0.83 | 0.08 | 4 |

![saq: inputs compared]({{ '/assets/img/sine-ordination/saq_inputs_v4.png' | relative_url }})

- With substitutions and indels the UMAP picture has three islands: {s1, s2, s7g, s9 and 53 s3 copies}, {s8 and 46 s3 copies}, {s4, s5, s6}.
- **Nested, substitutions only: 8 groups.** Published subfamily and its largest share in one group: s5 88 of 90, s8 95 of 100, s6 81 of 94, s1 82 of 98, s2 76 of 98, s4 72 of 88. One group holds all 99 s7g copies together with 65 s3 and 56 s9 copies. A further group of 35 copies is 32 s3. s9 has no group of its own; its copies sit with s7g (56), s1 (34) and s2 (9).
- **Nested, substitutions + indels: 12 groups**, among them s3 in two groups (51 of 65 and 42 of 47), s2 in two and s5 in two. ARI against the published labels 0.54 and 0.57 (K2P route: 0.19).
- The number of groups is not settled: 8 or 12 depending on indels, and the toy says the procedure tends to one or two too many. Several groups still contain co-varying characters that were not split (for example 25 in the group that is mostly s1).
- The label-free selection does as well as the columns where the published consensuses differ.

So the statement made here earlier, that no picture of copy distances can draw the saq subfamilies apart, was wrong: it held for K2P over all columns. On positions that co-vary, most published subfamilies are recovered without labels. Full tables: `saq_inputs_v4_crosstabs.txt`; every character with its weight, the consensus states and its frequency in each published subfamily: `saq_inputs_v4_features.tsv`.

Not done: the groups shown as alignments; the 54 unexplained columns; the cause of the over-splitting; a toy with recombination.

### 2.10 Four families, and the groups as alignments

Scripts: `projects/sine-ordination-toy/families/prep_family.py` and `family_inputs_v4.py`; method of 2.9 unchanged. Therioserver, 2026-11-12. Results per family in the same folder: `<name>_v4.json`, `_v4_crosstabs.txt`, `_v4_features.tsv`, `_v4.png`, `_prep.json`.

**Input.** The 100 random copies per subfamily from the Tal repository (`tim/alignments`, `rsi_final_20261102/alignments`; saq as in 2.8). Two *Timema cristinae* families and the six single-unit subfamilies of *Rhinolophus sinicus* that have 100 copies each (the composite subfamilies were left out: they are several units long and do not belong in one frame with single units). Body = the columns spanned by the subfamily consensus in the published alignment. MAFFT 7.526 `--auto` of consensuses + bodies.

**Flank gate, within each set** (method of 2.8): a copy was dropped if it shares both flanks above 60 % with another copy, or if a flank is shorter than 30 bp and cannot be tested. Dropped: *Timema* t6 127 of 500 (2 shared, 125 short), t1 78 of 400 (2 shared, 76 short), bat 1 of 600. The short flanks in *Timema* come from rows where the published alignment places part of the flank inside the consensus span. Not genome-wide in any set.

Each cell with three numbers: share of copies whose 15 nearest copies are mostly their own published subfamily / silhouette of the published labels / HDBSCAN groups on the UMAP picture.

| Family | Copies | Published subfamilies | K2P, CpG masked | Selected positions (substitutions + indels) | Consensus-difference columns selected | Split-half: selected in both halves, of per half | Nested groups, substitutions (ARI) | Nested groups, with indels (ARI) |
|---|---|---|---|---|---|---|---|---|
| *Scalopus aquaticus* Tal, s1–s9 | 866 | 9 | 0.64 / 0.01 / 2 | 0.83 / 0.07 / 3 | 80 of 87 | 78–99 of 89–139 | 8 (0.54) | 12 (0.57) |
| *Timema cristinae* t6-1 … t6-5 | 373 | 5 | 0.91 / 0.18 / 5 | 0.93 / 0.32 / 5 | 178 of 213 | 160–185 of 189–253 | 6 (0.71) | 5 (0.81) |
| *Timema cristinae* t1_1 … t1-4 | 322 | 4 | 0.71 / 0.15 / 3 | 0.81 / 0.31 / 4 | 30 of 72 | 7–22 of 23–63 | 7 (0.42) | 5 (0.41) |
| *Rhinolophus sinicus* units r5–r10 | 599 | 6 | 0.92 / 0.20 / 5 | 0.97 / 0.50 / 6 | 89 of 106 | 88–95 of 92–116 | 6 (0.88) | 9 (0.83) |

What the families show:

- **Bat units.** K2P already separates them fairly (5 groups). On selected positions the silhouette rises from 0.20 to 0.50 and the nested procedure (substitutions) returns 6 groups, one per published subfamily: r9 98 of 100, r10 100 of 100, r5 98, r6 82 of 99, r7 97, r8 92. With indels r6, r7 and r8 are each cut in two.
- ***Timema* t6.** K2P is already good here too. Nested with indels: 5 groups; t6-2, t6-3 and t6-4 each form one group (76 of 77, 55 of 61, 76 of 77); 22 of the 77 t6-1 copies sit with t6-5, the rest form their own group.
- ***Timema* t1.** The weakest case. t1_1 comes out as three groups of about 30 copies in both variants. t1-2 and t1-4 are not separated from each other. Only 30 of the 72 consensus-difference columns were selected and the split-half agreement is poor (7–22 shared). The set is small (322 copies, 161 per half), so the positions may be real and under-sampled, or not real; this run cannot tell.
- **saq** as in 2.9.
- Across the four, the selected-position distance raises the silhouette of the published labels in every family (0.01 → 0.07, 0.18 → 0.32, 0.15 → 0.31, 0.20 → 0.50) and the nearest-copy agreement in every family.

**The groups as alignments.** Rows: the published consensuses, then the copies ordered by nested group. Copy names: `g<group by substitutions>_i<group with indels>|<published label>|<coordinates>`. "Selected columns" shows only the positions that were selected or that separate two groups; "all columns" shows every column covered by at least half the copies. The groups are the procedure's candidates, not subfamily calls.

- *Scalopus aquaticus* Tal, s1–s9: [selected columns (137)](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2Fknowlege%2Fmain%2Fprojects%2Fsine-ordination-toy%2Ffamilies%2Fsaq_groups_selected.aln.fa&title=saq%20groups%2C%20selected%20columns) · [all columns (279)](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2Fknowlege%2Fmain%2Fprojects%2Fsine-ordination-toy%2Ffamilies%2Fsaq_groups_full.aln.fa&title=saq%20groups%2C%20all%20columns)
- *Timema cristinae* t6-1 … t6-5: [selected columns (244)](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2Fknowlege%2Fmain%2Fprojects%2Fsine-ordination-toy%2Ffamilies%2Ftim_t6_groups_selected.aln.fa&title=tim_t6%20groups%2C%20selected%20columns) · [all columns (427)](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2Fknowlege%2Fmain%2Fprojects%2Fsine-ordination-toy%2Ffamilies%2Ftim_t6_groups_full.aln.fa&title=tim_t6%20groups%2C%20all%20columns)
- *Timema cristinae* t1_1 … t1-4: [selected columns (103)](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2Fknowlege%2Fmain%2Fprojects%2Fsine-ordination-toy%2Ffamilies%2Ftim_t1_groups_selected.aln.fa&title=tim_t1%20groups%2C%20selected%20columns) · [all columns (269)](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2Fknowlege%2Fmain%2Fprojects%2Fsine-ordination-toy%2Ffamilies%2Ftim_t1_groups_full.aln.fa&title=tim_t1%20groups%2C%20all%20columns)
- *Rhinolophus sinicus* units r5–r10: [selected columns (102)](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2Fknowlege%2Fmain%2Fprojects%2Fsine-ordination-toy%2Ffamilies%2Frsi_units_groups_selected.aln.fa&title=rsi_units%20groups%2C%20selected%20columns) · [all columns (250)](https://toki-bio.github.io/MSA-viewer/?url=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2Fknowlege%2Fmain%2Fprojects%2Fsine-ordination-toy%2Ffamilies%2Frsi_units_groups_full.aln.fa&title=rsi_units%20groups%2C%20all%20columns)

![tim_t6: inputs compared]({{ '/assets/img/sine-ordination/tim_t6_v4.png' | relative_url }})

![tim_t1: inputs compared]({{ '/assets/img/sine-ordination/tim_t1_v4.png' | relative_url }})

![rsi_units: inputs compared]({{ '/assets/img/sine-ordination/rsi_units_v4.png' | relative_url }})

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

The ordination rows (first five) were done on 2026-11-11 in `ordination_toy_v2.py`; results in [2.6](#26-repaired-run-5-seeds). The last two rows (pan-SINEome toy) are still open; the Jaccard negative-eigenvalue share was only recomputed once and found to be 0.

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

Indels, 5′ truncation and CpG decay: done in [2.7](#27-harder-toy-indels-truncation-cpg-decay-number-of-subfamilies-not-given). ARI fell below 0.99 as required, and the methods now separate: K2P with pairwise deletion + t-SNE or UMAP holds; one-hot does not.

Still to add, each being a named cause of label disagreement on real data (Carey et al. 2021): recombination / gene conversion between subfamilies, a continuum of copy ages without discrete groups, and very unequal subfamily sizes (a small subfamily next to a large one).

### 4.3 One real family from a SINEderella run

- Input: the copies of one family that the workflow accepted (unique flanks), in the alignment the workflow already makes. Nothing new is searched or aligned.
- Distance (supported by the toy in 2.7): K2P with pairwise deletion over shared columns, CpG columns of the family consensus masked. Not one-hot. Also exclude the poly-A tail and low-complexity 3′ parts, as the Alu network studies and the bat Ves SINE study did (reported by the literature agent from full texts; not in section 5).
- Picture: t-SNE or UMAP of that distance, with PCoA beside it, colored by the subfamily the workflow assigned.
- Question: do the assigned subfamilies appear as groups, and which copies sit between groups?
- Control: a family whose subfamilies are already established should show them. If it does not, the picture is not informative for this kind of data and the idea stops here.
- Possible check without known labels (after Carey et al. 2021, who used segmental duplicates and human–chimp orthologs): orthologous copies of one insertion in two related genomes, from SINE_orth_loc, should get the same subfamily. Segmental duplicates are excluded here, since their flanks are shared.

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

**Subfamily delimitation in practice, and simulation (sources for 2.7)**

- Zhou Y, et al. (2019) Subtype classification and functional annotation of L1Md retrotransposon promoters. *Mobile DNA* 10:14. [doi:10.1186/s13100-019-0156-5](https://doi.org/10.1186/s13100-019-0156-5). Full text. Truncated monomers excluded ("we defined a set of 'core' monomers ... Monomers not meeting these criteria were excluded"); PCA to 100 axes, HDBSCAN (min points 10–20), t-SNE for display only.
- Storer J, Hubley R, Rosen J, Smit AF (2021) Curation guidelines for de novo generated transposable element families. *Curr Protoc* 1:e154. [doi:10.1002/cpz1.154](https://doi.org/10.1002/cpz1.154). Full text. COSEG on full-length copies; "aligned TG and CA dinucleotides sequences are the result of a single differential decay of a CpG site and do not constitute co-segregating individual substitutions"; "COSEG is not appropriate because depending upon the number of sequences that appear 'truncated', a large section of sequences will be missed."
- Carey KM, Patterson G, Wheeler TJ (2021) Transposable element subfamily annotation has a reproducibility problem. *Mobile DNA* 12:4. [doi:10.1186/s13100-021-00232-4](https://doi.org/10.1186/s13100-021-00232-4). Full text. Two copies of one ancestral insertion get different subfamily labels in over 10 % of cases (Alu 12.4 %, L1 14.1 %, MIR 10.4 %); recombination and chance mutation at diagnostic sites are the causes discussed.
- Wacholder AC, et al. (2014) Inference of transposable element ancestry. *PLoS Genet* 10:e1004482. [doi:10.1371/journal.pgen.1004482](https://doi.org/10.1371/journal.pgen.1004482). Full text; author list not checked.
- Vassetzky NS, Kosushkin SA, Korchagin VI, et al. (2021) New Ther1-derived SINE Squam3 in scaled reptiles. *Mobile DNA* 12. [doi:10.1186/s13100-021-00238-y](https://doi.org/10.1186/s13100-021-00238-y). Full text. Subfamilies told apart by insertions/deletions, because "single-nucleotide mutations can be used to identify subfamilies only for highly conserved SINEs".
- Hubley R, Wheeler TJ, Smit AF (2022) Accuracy of multiple sequence alignment methods in the reconstruction of transposable element families. *NAR Genom Bioinform* 4:lqac040. [doi:10.1093/nargab/lqac040](https://doi.org/10.1093/nargab/lqac040). Full text. Simulated TE copies: indel rate 0.20 of substitutions (insertion 0.08, deletion 0.12), length mean 1.7, max 20. SINE simulation only planned: "additionally consider SINE and LTR TE classes".
- Miao, Li (2026) TEvarSim (initials not checked). *PLoS Comput Biol* 22:e1013933. [doi:10.1371/journal.pcbi.1013933](https://doi.org/10.1371/journal.pcbi.1013933). Full text. Defaults: proportion of sequences truncated 0.3, maximum proportion truncated 0.5.
- Xing J, et al. (2004) *J Mol Biol* 344:675-682. [doi:10.1016/j.jmb.2004.09.058](https://doi.org/10.1016/j.jmb.2004.09.058). Abstract: "a relatively constant CpG versus non-CpG substitution ratio of approximately 6 for the young (AluY) and intermediate (AluS) Alu subfamilies".
- Lunter G, Ponting CP, Hein J (2006) *PLoS Comput Biol* 2:e5. [doi:10.1371/journal.pcbi.0020005](https://doi.org/10.1371/journal.pcbi.0020005). Indels "approximately 10-fold less frequent" than substitutions (quote as returned by the fetch tool; full text not read).

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

Run the scripts in `projects/sine-ordination-toy/`; the first two write their figures to the working directory. Versions used for the figures and the first table: Python 3.13, numpy, scipy, scikit-learn, matplotlib, umap-learn 0.5.12, pacmap 0.9.1. `ordination_toy_v2.py [NSEED]` (versions in 2.6) writes `ordination_v2_metrics.json` and no figures. `ordination_toy.py` also writes `ordination_metrics.json`; `pansineome_toy.py` writes `toy_pansineome_PAV.tsv`.
