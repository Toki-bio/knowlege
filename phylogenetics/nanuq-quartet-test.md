---
title: NANUQ quartet test for newbies
type: explanation
---

> Type: explanation

# NANUQ: the quartet test, without the usual muddle

NANUQ is a method for asking whether multilocus gene trees look like they came from a **species tree** plus incomplete lineage sorting, or from a **species network** with hybridization / gene flow.

Paper: Allman, Baños & Rhodes (2019), *NANUQ: a method for inferring species networks from gene trees under the coalescent model*. Implemented in the R package [MSCquartets](https://cran.r-project.org/package=MSCquartets).

The acronym is **N**etwork inference **A**lgorithm via **N**eighbourNet **U**sing **Q**uartet distance.

One sentence:

> For every set of four taxa, count which of the three unrooted gene-tree shapes appears, test whether those counts fit a tree, a star, or a 4-cycle, then turn the tests into a distance that NeighborNet can draw as a splits graph.

It does **not** take raw sequences as input. It takes a **collection of gene trees**.

---

## Why four taxa, not three

Three taxa **can** show different branching. The question is *which kind of tree you are looking at*.

### Rooted trees (time has a direction)

A **root** is a common ancestor that orients time: who split first.

On three taxa there are three different rooted trees (rooted triples):

```text
(A,(B,C))     (B,(A,C))     (C,(A,B))
```

Those are genuinely different histories. Methods that use rooted triples are testing something real. So “three is not enough for branching” is false.

### Unrooted trees (no start point)

Most gene-tree summaries used by NANUQ are **unrooted**. There is no marked ancestor. You only have how the living taxa connect.

An unrooted binary tree on *n* taxa has **n − 3 internal branches**:

| Taxa | Internal branches | Unrooted shapes |
| --- | --- | --- |
| 2 | −1 (just a line) | 1 |
| 3 | 0 | 1 (a Y / three-way junction) |
| 4 | 1 | 3 competing quartets |

On three unrooted taxa there is only a single junction. Cut any branch and you isolate one taxon. There is no internal branch, and no second unrooted shape to reject.

On four taxa you finally get one internal branch, and three mutually exclusive unrooted topologies:

```text
AB | CD          AC | BD          AD | BC
```

That internal branch is the thing a quartet test can measure. NANUQ lives in this unrooted world, so **four is the minimum for this test**, not the minimum for evolution.

### What “start point” means

The start point is the **root**: the ancestor that gives time a direction.

Without it, the rooted histories `(A,(B,C))` and `((A,B),C)` collapse to the same unrooted Y. You can still see that some pairs are genetically closer (branch lengths / distances). You cannot, from unrooted three-taxon geometry alone, choose between competing unrooted topologies, because there is only one.

Similarity is not the same as a testable unrooted split. Relatedness is defined by shared ancestry (timing of splits), not by who still has wings.

---

## What NANUQ actually does

### Input

Many inferred **gene trees** (one per locus), possibly with missing taxa. Each subset of four taxa must appear on at least one gene tree.

### Step 1 — count quartets

For taxa A, B, C, D, look at every gene tree that contains them and tally the three resolved shapes:

```text
qcCF = (count of AB|CD,  count of AC|BD,  count of AD|BC)
```

Normalize the counts to frequencies and you have an **empirical quartet concordance factor** (CF). It is a point in a triangle (a 2-simplex) whose three corners are “100% of gene trees show this one shape.”

Under a species **tree** plus incomplete lineage sorting (the multispecies coalescent, MSC):

* one topology is most frequent (the species-tree quartet);
* the other two are **equal** (the two mismatch topologies produced by ILS).

That is the “tree-like” signature: one big count, two smaller equal counts. See also [coalescence as same chromosome copy](https://toki-bio.github.io/knowlege/population-genetics/coalescence-same-chromosome-copy/).

If the four taxa sit on a **4-cycle** (hybridization around those four), the two minor counts are **not** equal in the way the tree MSC requires. The CF point is pulled off the tree-model lines toward an edge of the triangle.

If all three counts are about 1/3, the quartet is consistent with a **star** (no internal branch; a polytomy), or with ILS so extreme that the internal branch is effectively zero.

### Step 2 — hypothesis tests

NANUQ does not eyeball the triangle and stop. For each quartet it tests, roughly:

* tree (MSC on a resolved quartet);
* star (unresolved);
* 4-cycle / hybridization (reject the tree MSC in the direction expected for a cycle).

Two significance levels matter in MSCquartets: **α** (tree vs not-tree / hybridization) and **β** (star). Chan et al. (below) used a very small α (10⁻⁷⁰) so that calling a 4-cycle was conservative. They still found many 4-cycle quartets.

### Step 3 — distance, then a picture

Test results are converted into a **NANUQ quartet distance** between taxa. NeighborNet (often in SplitsTree) draws a **splits graph**.

* Tree-like data: mostly ordinary branches.
* Reticulation: boxes / grids / extra parallel edges.

The splits graph is not automatically “the true network.” NANUQ is statistically consistent for **level-1** networks under the network MSC (at most one cycle per blob). Human interpretation of the tests **and** the graph is part of the method. If the data are not a level-1 network, a pretty NeighborNet can still be drawn; that does not make it true.

SNaQ is a related CF-based network method. NANUQ is often used to see *whether* quartets are tree-like before or beside a full network search.

---

## How to read the triangle (simplex plot)

Each **dot is one set of four taxa**, not one gene.

### Geometry (where the point sits)

```text
            AB|CD
           /     \
          /       \
         /  center \     center ≈ (1/3, 1/3, 1/3)  star / extreme ILS
        /           \
   AC|BD ----------- AD|BC
```

* **Near a corner:** almost all gene trees show one quartet. Clean, tree-like, long internal branch.
* **On a model line from the center toward a corner:** tree MSC: major topology plus two equal minor ones (ILS).
* **Off those lines, toward an edge:** two topologies compete in a way a tree MSC does not predict — the hybridization / 4-cycle signal NANUQ is after.
* **At the center:** all three topologies equally common.

An edge of the triangle means only two of the three quartets appear. That is *not* the same as “a triangle side = a terminal branch to species C.” Branches of the species tree are not the sides of this plot.

### Chan et al. symbol coding (Fig. 3a, 5a)

In [Chan et al. 2026](https://doi.org/10.1093/sysbio/syag001) on Bornean fanged frogs (*Limnonectes*):

* **circles** = quartets accepted as a bifurcating tree;
* **triangles** = quartets classified as a 4-cycle (their “gene flow” call).

So you read two layers: *where* the point is, and *what the test decided*.

---

## NANUQ vs ASTRAL vs ABBA–BABA

All three can use four taxa. They answer different questions.

| | **ABBA–BABA (D)** | **ASTRAL** | **NANUQ** |
| --- | --- | --- | --- |
| Input | SNPs / site patterns | gene trees | gene trees |
| Needs a root / outgroup? | Yes: `(((P1,P2),P3),O)` | No | No |
| ILS | Background: ABBA ≈ BABA under ILS only | Allowed (MSC on a **tree**) | Allowed (MSC / network MSC) |
| Hybridization | The thing being tested (excess ABBA or BABA) | Not the target; forced into a tree | The thing being tested (4-cycle CFs) |
| Output | A statistic (and p-value) for named taxa | One species **tree** | Tests + distance → **splits graph / network** |

ASTRAL searches for the species tree that matches the largest number of gene-tree quartets. Disagreement among genes is treated as ILS around a tree, not as evidence for extra edges.

ABBA–BABA does not build a tree or a network. It tests whether P3 shares more derived alleles with P2 than with P1, given an outgroup.

NANUQ does not replace either. It asks, quartet by quartet, whether the CF vector is tree-like under the coalescent or looks like a cycle, then assembles those decisions.

---

## Worked example: Chan et al. fanged frogs

Chan, Neokleous, Anuar, Brown, Hutter, Das & Hertwig (2026), *Systematic Biology*. [doi:10.1093/sysbio/syag001](https://doi.org/10.1093/sysbio/syag001)

Mitochondrial trees had split Bornean *Limnonectes kuhlii*-complex frogs into many cryptic “species.” Nuclear target-capture (~13k loci) plus methods that allow gene flow collapsed that to **six or seven independently evolving lineages**. Genome-wide net divergence sat in a **0.5–2% gray zone** of the speciation continuum.

Useful terms from the paper:

* **Artifactual branch effect** — gene flow can make an admixed group look like a long, early-diverging branch on a tree, as if it were ancient and isolated.
* **Species-definition anomaly zone** — distances within a “species” can exceed distances between “species” when gene flow is still happening.

### Which markers contradicted UCE?

Not a vague “SNPs vs UCEs” slogan. They ran NANUQ **separately** on marker classes (their Fig. 5a):

* **UCE:** no reticulation in the NANUQ simplex (all tree-like quartets under their threshold).
* **AHE, BUSCO, FrogCap:** variable numbers of 4-cycle quartets.

On recent migration (BPP MSC+M, their Fig. 5b):

* **AHE and UCE:** no significant recent migration (M < 0.1).
* **BUSCO and FrogCap:** migration detected for some population pairs.

Dsuite / f-branch also found introgression, including older events that UCE/AHE can still see. The authors’ own summary: AHE and UCE are useful for deep history and **older** introgression, and are **poor at recent gene flow**. That is why a UCE-only tree can look clean while other markers say the lineages still mix.

Population-structure SNPs (sNMF, t-SNE) were a different analysis: they asked which samples cluster, not which quartets are 4-cycles.

---

## Would NANUQ help SINE subfamily discrimination?

Short answer: the **quartet idea** is useful; dropping NANUQ into [SINEplot.py](https://github.com/Toki-bio/SINEplot/blob/main/SINEplot.py) as-is is not.

SINEplot places copies from **bitscores** (ssearch) with MDS. MDS is pairwise similarity. Pairwise plots cannot tell a clean split from a recombinant / intermediate copy: both can sit “between” two clouds.

NANUQ expects **orthologous gene trees** generated under a coalescent on a species tree or network. SINE copies are mostly **paralogs** from master-copy retrotransposition. Treating every copy as a “taxon” and every pairwise alignment as a “gene tree” is a different statistical object. A bitscore ratio of 2nd-best / 1st-best is an ambiguity heuristic. It is **not** a concordance factor and **not** NANUQ.

### Combinatorial cost

Number of quartets is C(*n*, 4). About 4 million for *n* = 100; impossible for typical copy catalogs. Full NANUQ on all copies is the wrong default.

### What is actually implementable

1. **Subfamily consensi as taxa.** If you have 6–12 named subfamilies, C(12, 4) = 495. Build one tree (or a few locus-like partitions) on those consensi and *then* a quartet / network method is tractable. That tests whether **named subfamilies** form a tree or a recombination network, not whether every copy is pure.
2. **Random or targeted quartets**, not all of them. Sample; or only test copies whose top two subfamily bitscores are close, plus two consensi and a distant outgroup consensus.
3. **Do not color MDS dots and call it NANUQ.** If you flag unstable copies, say so: “top-two bitscores too close,” not “4-cycle under the network MSC.”

The Chan lesson that *does* transfer: pairwise divergence and a forced tree will over-split. Intermediate SINEs and “gray zone” subfamilies need a conflict diagnostic, not a new cluster label for every 1% of sequence difference.

Related notes in this repo: [clustering for SINE subfamily discovery](https://toki-bio.github.io/knowlege/methods/clustering-for-sine-subfamily-discovery/), [PCA vs UMAP/t-SNE](https://toki-bio.github.io/knowlege/methods/dim-reduction-for-sine-separation/).

---

## Bottom line

```text
3 taxa + a root  →  real competing histories (rooted triples)
3 taxa, unrooted →  one Y; nothing for NANUQ to test
4 taxa, unrooted →  three quartets and one internal branch

NANUQ counts those quartets across gene trees,
tests tree vs star vs 4-cycle,
draws a splits graph.

ASTRAL: same quartets, forced into a tree.
ABBA–BABA: rooted four-taxon SNP test for gene flow.
SINE copies: use quartet thinking; do not pretend bitscore MDS is NANUQ.
```
