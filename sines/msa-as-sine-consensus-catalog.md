---
title: MSA as the SINE consensus catalog
type: explanation
---

> Type: explanation

# The MSA is the SINE consensus catalog

## The problem this note is really about

You have accumulated many **SINE consensus sequences from many taxa**.

They look roughly like:

```text
taxon A ── SINE family X ── several related consensuses
taxon B ── SINE family Y ── somewhat related to X
taxon C ── SINE family Z ── distant
...
```

Some consensuses are obviously close. Others are very distant. Between those extremes sit the hard cases:

> Same family?  
> Different subfamilies?  
> Distantly related families?  
> Or unrelated SINEs?

And the point is not to have a clustering algorithm announce the answer. The point is to **keep enough evidence that you can define and revise those boundaries yourself**.

So the computational problem is not primarily:

> How do I invent a clever mathematical object to replace the alignment?

It is:

> How do I maintain a collection of complete SINE sequences and expose their relationships well enough that biological group boundaries can be drawn, inspected, and revised?

That distinction decides the architecture.

---

## The wrong turn: inventing a structure to replace the MSA

A natural temptation, especially if one has pangenome graphs in mind, is to ask:

> What graph / POA / compact structure should store a SINE family?

That question is seductive and usually premature.

A variation graph is a representation *derived from* relationships among sequences. It can be useful later. It is not required for the fundamental storage problem, and it risks "mincing" the thing you care about most: each consensus as a complete sequence path.

What you need first is much less exotic.

---

## The key insight: every consensus is already a path through the alignment

Suppose you have:

```text
                  111111111122222222223333333333
                  123456789012345678901234567890
SINE_001          ATGCGT---AATCCGGT---TACCG...
SINE_002          ATGCGT---AATTCGGT---TACCG...
SINE_003          ATG-GTCC-AATTCGGT---TACCG...
SINE_004          GGGACCAA-CCGTAGC----GTTAA...
```

The scientifically important information is not only column summaries such as:

> position 17 = 83% A

It is also:

> **SINE_002 is one complete sequence that differs from SINE_001 here and resembles SINE_003 there.**

The MSA already preserves that.

Every row is a path through a shared coordinate system. Substitutions, insertions, deletions, and the co-occurrence of changes within one sequence are all still there. You can open the file in an ordinary alignment viewer. No specialized graph infrastructure is required.

### Therefore:

**The MSA itself is your lightweight pangenome representation.**

The "pangenome-like" property comes from retaining **multiple complete sequence paths through a shared coordinate system**, not from constructing a formal variation graph.

---

## Make the central object an alignment collection

Rather than a "SINE graph," the central object should be something like a **SINE Alignment Collection**.

Conceptually:

```text
collection/
    sequences.fa
    alignment.fa
    metadata.tsv
    distances.tsv
    tree.nwk
    groups.tsv
    parameters.yaml
```

### `sequences.fa` — primary archive

The untouched consensuses:

```text
>SINE000183 taxon=Mus_musculus ...
ATGCG...
```

This is the source of truth for sequence content. Nothing derived should replace it.

### `alignment.fa` — scientific core

The MSA containing exactly those sequences.

This is probably the most important object scientifically because it preserves:

* substitutions,
* insertions and deletions,
* relative positions,
* complete sequence identity,
* relationships among changes occurring in the same sequence.

### `metadata.tsv` — revisable annotations

```text
id       taxon       source       family    subfamily    group
SINE001  Homo        RepBase      F001      F001.1       ...
SINE002  Homo        custom       F001      F001.1       ...
SINE003  Macaca      custom       F001      F001.2       ...
```

Family and subfamily labels are **annotations on sequences**, not properties baked into the representation. That matters because you will revise them.

---

## Separate measurement from interpretation

From the MSA, compute a pairwise relationship matrix:

|       | SINE1 | SINE2 | SINE3 | SINE4 |
| ----- | ----: | ----: | ----: | ----: |
| SINE1 |     0 |   .04 |   .07 |   .31 |
| SINE2 |   .04 |     0 |   .05 |   .29 |
| SINE3 |   .07 |   .05 |     0 |   .27 |
| SINE4 |   .31 |   .29 |   .27 |     0 |

Then keep two layers separate:

### Measurement

```text
distance(SINE_A, SINE_B) = 0.073
```

### Interpretation

```text
0.073 → same subfamily   (under threshold T, method D, dataset version V)
```

The first is data.  
The second is your biological classification.

That separation makes the catalog revisable and scientifically defensible. Two years later, when 2,000 new consensuses arrive, you can still say how `F001.1` was defined.

---

## The colored matrix is the decision interface

This is where the "colored alignment / shaded relatedness" intuition becomes practical — not as a replacement for the MSA, but as a view of whole-sequence relationships derived from it:

```text
             sequence
          1  2  3  4  5
       1  █  ▓  ▓  ░  ░
       2  ▓  █  ▓  ░  ░
       3  ▓  ▓  █  ░  ░
       4  ░  ░  ░  █  ▒
       5  ░  ░  ░  ▒  █
```

Blocks emerge visually. Selecting a block can mean:

> These sequences look like a group. Assign `SUBFAMILY = F001.2`.

The alignment and original sequences never disappear.

A useful interface sketch:

```text
Upper panel:   similarity / distance heatmap
Side panel:    current family / subfamily tree of labels
Lower panel:   the actual MSA for the selected sequences
```

That combination is much closer to the real workflow than either "only a tree cut" or "only a graph."

---

## Do not let clustering become the authority

A dendrogram is useful. Automatic cuts are dangerous as definitions.

Biological relationships are not guaranteed to be cleanly hierarchical. Consider:

```text
A ─ B ─ C ─ D
```

where A–B, B–C, and C–D are each moderate, but A–D is large. Single-linkage clustering can then invent absurd "families" by **chaining**.

So the system should expose, side by side:

1. pairwise distances,
2. the MSA,
3. trees / clustering as optional views,
4. manually assigned groups.

Clustering may propose. It should not define. That is the same principle as in the [clustering for SINE subfamily discovery](https://toki-bio.github.io/knowlege/methods/clustering-for-sine-subfamily-discovery/) note: unsupervised methods are hypothesis generators, not taxonomies.

Thresholds such as:

```text
below X = same subfamily
below Y = same family
below Z = same broader lineage
```

are your classification layer. They are not intrinsic properties of the storage format.

---

## Distance must keep indel information — but that does not require a graph

"Preserve indel information" does not mean "build a graph." It means the **distance function must not throw indels away**.

Candidates to compare against expert judgment:

### Simple alignment distance

In homologous columns:

```text
A vs A = 0
A vs G = 1
A vs - = 1
- vs - = ignored
```

Transparent and easy to audit.

### Gap-event distance

Treat a run of gaps as one indel event rather than many independent column differences. Often closer to biological intuition.

### Hybrid

```text
D = substitution component + weighted indel component
```

Do not choose this theoretically. Benchmark several definitions against sequences whose relationships you already understand:

```text
same subfamily
same family, different subfamily
related family
unrelated
```

If a distance produces heatmaps whose blocks match your judgment, that distance becomes the quantitative backbone of the catalog. If it does not, adjust the distance — do not invent another storage architecture.

---

## Alignment quality forces a two-level system

The distance matrix is only as meaningful as the MSA.

If two very divergent SINEs are forced into one bad alignment, the matrix will quantify the bad alignment with great confidence.

So distinguish:

### Within a reasonably homologous group

Use a shared MSA and MSA-derived distances.

### Across very divergent SINEs

Do not force everything into one giant alignment. First decide whether sequences belong in the same comparison space at all.

A realistic catalog hierarchy is therefore operational rather than metaphysical:

```text
GLOBAL SINE CATALOG
        │
        ├── search / relationship space
        │
        ├── broad homologous groups
        │       ├── alignment
        │       ├── distance matrix
        │       ├── tree
        │       └── family / subfamily assignments
        │
        └── other groups ...
```

A distant SINE does not need to participate in the same MSA as everything else.

---

## Put HMMs in the search layer, not the storage layer

The intuition that an HMM is insufficient as the **storage representation** is correct. An HMM summarizes a family; it does not retain every complete consensus as a first-class object.

But HMMs are excellent as a **search representation**.

```text
                 SOURCE OF TRUTH
              sequences + MSA
                     │
          ┌──────────┼──────────┐
          │          │          │
       distance     tree       HMM
          │                     │
    classification           genome scan
```

The HMM is a projection. If search later confuses `F001` and `F002`, you revise the search model without rewriting the underlying catalog.

---

## Multiple representations of one collection

You do not need a mysterious intermediate between "pangenome" and "HMM."

You need several representations of the same sequence collection, each with a job:

| Representation | Purpose |
| --- | --- |
| FASTA | original consensuses |
| MSA | exact positional relationships + indels |
| distance matrix | whole-sequence relatedness |
| dendrogram / tree | relationship visualization |
| family / subfamily labels | biological interpretation |
| HMM | genome search |
| representative sequences / medoids | fast search and classification |
| optional graph later | specialized analysis, if needed |

The MSA is the connective tissue. Everything else hangs off it.

---

## Useful extras that still do not replace the members

### Group representatives

After selecting a group, you can derive:

* a **medoid** — the actual member with smallest average distance,
* an ordinary consensus,
* several representatives,
* an ancestral reconstruction, if useful.

None of these replaces the members. They are conveniences for search and display.

### Nearest-neighbor signatures

Store for each consensus:

```text
SINE001
  nearest:
    SINE002   0.031
    SINE014   0.037
    SINE008   0.044
```

Classifying a new consensus then asks:

> Which existing SINEs is this actually closest to?

rather than only:

> Which HMM does this hit?

### Relationship atlas

Dimensionality reduction (UMAP / t-SNE) can provide a global visual atlas of consensuses. Color by assigned family, shape by taxon, border by subfamily. Selecting a point opens the sequence, its alignment neighborhood, nearest neighbors, and distance to current boundaries.

Important constraint:

**UMAP does not define the taxonomy.**  
It visualizes relationships already measured in the distance matrix.

---

## What to build first

At this stage, do not build a graph database.

Build a boring, robust catalog:

```text
SINEdb/
├── sequences.fa
├── metadata.tsv
├── alignments/
│   ├── group_001.fa
│   └── ...
├── distances/
│   ├── group_001.tsv
│   └── ...
├── trees/
│   └── ...
├── classifications.tsv
├── thresholds.yaml
└── visualizations/
```

With classifications like:

```text
sequence_id   broad_group   family   subfamily   status
SINE001       G001          F001     F001.1      accepted
SINE003       G001          F001     F001.2      accepted
```

None of these labels modifies the sequence or the alignment. That is the whole point.

---

## Core diagram

```text
                     SINE CONSENSUSES
                            │
                            ▼
                    ┌───────────────┐
                    │      MSA      │
                    │ complete seqs │
                    │ + indels      │
                    │ + connectivity│
                    └───────┬───────┘
                            │
                ┌───────────┴───────────┐
                ▼                       ▼
       PAIRWISE DISTANCES             TREE
                │                       │
                └───────────┬───────────┘
                            ▼
                     COLORED MATRIX
                            │
                            ▼
                  YOUR BIOLOGICAL CUTS
                            │
                 ┌──────────┼──────────┐
                 ▼          ▼          ▼
              FAMILY    SUBFAMILY   GROUP
                 │
                 ▼
        search projections:
        HMM / representatives / k-mers
```

* **MSA** = storage and relationship core  
* **colored matrix** = human decision interface  
* **classifications** = annotations  
* **HMMs** = search projections  

That is the architecture that matches the biological problem.

---

## The first empirical step before implementation

Take 50–200 consensuses whose relationships you already understand and mark:

```text
same subfamily
same family, different subfamily
related family
unrelated
```

Then test several indel-aware distances:

```text
your classification
        │
        ▼
 distance matrix → heatmap → threshold
        │
        ▼
 does it reproduce your judgment?
```

If yes, you have the quantitative backbone of the catalog.  
If no, change the distance definition — not the storage philosophy.

---

## Key takeaways

* The real problem is maintaining complete SINE consensuses and exposing relationships for **human-defined** family boundaries.
* Do not replace the alignment with a graph at the start. The MSA is already a lightweight pangenome of complete sequence paths.
* Keep sequences, MSA, distances, trees, labels, and HMMs as separate layers with separate jobs.
* Separate measurement (`distance = 0.073`) from interpretation (`same subfamily`).
* Use the colored distance matrix as the decision interface; keep the MSA visible underneath.
* Do not let dendrogram cuts become the definition of families.
* Preserve indels in the distance function; choose the function empirically.
* Use separate alignments for broad homologous groups rather than one forced global MSA.
* HMMs are for genome search, not for being the source of truth.
