---
title: HDBSCAN vs vsearch for SINE grouping
type: explanation
---

> Type: explanation

# HDBSCAN vs vsearch for SINE lineage grouping

## Short answer

**Yes — HDBSCAN can be a good choice for SINE grouping**, often better than treating t-SNE or PCA alone as a clustering method. But only in the right pipeline, and not as a replacement for sequence-identity clustering.

The useful comparison is:

| Tool | Question it answers |
| --- | --- |
| **vsearch** | Are these sequences ≥ X% identical? |
| **HDBSCAN** | Do these sequences form natural density groups in feature space? |

They solve **different layers** of the problem. The strongest practical answer is usually to use both.

This note sits next to the earlier discussion of [K-means / GMM as candidate generators](https://toki-bio.github.io/knowlege/methods/clustering-for-sine-subfamily-discovery/) and the [MSA-centered catalog](https://toki-bio.github.io/knowlege/sines/msa-as-sine-consensus-catalog/) idea: clustering may propose groups; it should not silently become the taxonomy.

---

## Why HDBSCAN fits SINE-shaped data

SINE families often show:

* uneven densities — dense young cores, diffuse old clouds;
* noise / borderline copies — truncated, chimeric, or badly assigned elements;
* non-linear structure in feature space;
* gradual subfamily divergence rather than clean spherical blobs.

HDBSCAN is attractive because it:

* does **not** require predefining the number of clusters;
* handles **variable-density** clusters;
* explicitly labels **noise**, which is valuable for junk or ambiguous SINEs.

That is a real advantage over:

* **PCA alone** — useful linear preprocessing / visualization, not a clustering method;
* **t-SNE** — useful for plots, unreliable as a clustering substrate because it distorts global distances and can invent compact blobs.

---

## Recommended pipeline

Do **not** run HDBSCAN on raw sequences as strings.

### 1. Feature space

Represent each SINE copy / consensus numerically, for example:

* k-mer frequencies (often 4–6-mers),
* alignment-derived similarity / embedding features,
* diagnostic-mutation or divergence features,
* existing pipeline metrics (e.g. SINEderella-style assignment features), if they encode biology rather than run artifacts.

Feature choice decides what biology you are allowed to discover. Bad features produce confident nonsense.

### 2. Optional reduction

Before HDBSCAN:

* light **PCA** (e.g. 20–50 dimensions), or
* **UMAP**, which is usually preferable to t-SNE if a non-linear embedding is needed for clustering.

Avoid:

```text
t-SNE → HDBSCAN
```

t-SNE is for looking. It is a poor coordinate system for defining families.

### 3. Clustering

Run HDBSCAN in that reduced feature space.

Key knobs:

* `min_cluster_size` — family / subfamily granularity;
* `min_samples` — how easily points become noise versus cluster members.

These are biological tuning parameters disguised as software settings. Document them.

---

## Why not t-SNE + clustering?

t-SNE:

* destroys global distances,
* can create artificial clusters,
* is highly sensitive to perplexity,
* is not designed as input to a density clusterer that you then treat as taxonomy.

Use t-SNE for figures. Use PCA / UMAP + HDBSCAN, or no embedding at all if the feature space is already compact and meaningful.

---

## When HDBSCAN fails

HDBSCAN is not magic.

It fails or misleads when:

* the feature space is poor;
* families are extremely overlapping (especially very old degraded SINEs);
* the dataset is too small for density estimation;
* you interpret noise labels as biological truth without looking at the sequences;
* you treat cluster IDs as final subfamily names.

As with GMM / K-means: **cluster → validate**, never **cluster → publish name**.

---

## HDBSCAN vs vsearch

### vsearch: sequence-identity clustering

Typical use:

```text
vsearch --cluster_fast / --cluster_size
--id 0.80–0.95
```

**Pros**

* deterministic and reproducible;
* biologically explicit percent-identity threshold;
* very fast and scalable;
* standard in TE / SINE workflows;
* excellent for consensus building and library compression.

**Cons**

* hard threshold forces an arbitrary cut;
* struggles with fragmented copies, length variation, and subfamily gradients;
* high identity over-splits;
* low identity over-merges.

SINE reality is often continuous divergence. A single identity cutoff is a blunt instrument.

### HDBSCAN: feature-space density clustering

**Pros**

* no fixed percent-identity cutoff;
* can reveal dense young cores vs diffuse old clouds;
* can separate hidden substructure;
* noise labeling helps quarantine ambiguous copies.

**Cons**

* depends on feature engineering;
* clusters are less explicitly “% identity = X”;
* harder to debug;
* less standard in TE pipelines, so validation burden is higher.

### Intuition

Imagine a SINE lineage with gradual mutations:

```text
VSEARCH at 90%
   → arbitrary hard splits along a continuum

HDBSCAN in a good feature space
   → dense young core
   → diffuse older cloud
   → maybe 2–3 real density groups + noise
```

Neither output is automatically the biological truth. Both are proposals.

---

## What works best in practice: hybrid roles

### Option A — classical and robust

```text
vsearch (e.g. 80–85%)
   → rough families
   → refine by alignment / consensus / diagnostics
```

Safe and publishable. May miss subtle substructure.

### Option B — feature-space structure

```text
features → PCA/UMAP → HDBSCAN
```

Richer structure, more validation needed.

### Option C — recommended combined use

Give each tool a job:

```text
VSEARCH
  → initial grouping
  → consensus building
  → redundancy reduction

HDBSCAN
  → look inside / across groups for
       hidden subfamilies
       misassigned elements
       borderline / noise cases
```

Then force every interesting HDBSCAN proposal back through:

* MSA / diagnostics,
* consensus recoverability,
* phylogenetic coherence,
* and, if relevant, the [MSA catalog logic](https://toki-bio.github.io/knowlege/sines/msa-as-sine-consensus-catalog/) of measurement versus interpretation.

---

## Placement relative to SINEderella-style outputs

If a pipeline already produces assignment tables and locus beds, HDBSCAN is usually most useful as a **second-pass audit**, not as the primary assigner.

Example roles:

* take copies assigned to one broad family and ask whether feature space contains multiple dense cores;
* take “noise” or low-confidence points and inspect whether they are junk, chimeras, or real rare lineages;
* compare HDBSCAN groups against existing canonical labels to find disagreements worth manual review.

Do not silently overwrite a curated assignment system with HDBSCAN cluster IDs. Especially remember that different SINE classification systems are not interchangeable; unsupervised clusters are yet another system unless validated into the one you mean.

---

## Bottom line

* **PCA** — good preprocessing / exploration.
* **t-SNE** — visualization only.
* **vsearch** — standard, robust, percent-identity grouping and consensus workflows.
* **HDBSCAN** — strong density-based proposal method for uneven SINE structure and noise, after sensible features.
* **Best practical answer** — vsearch for initial grouping / consensus; HDBSCAN for discovering hidden structure and suspicious cases; biology (alignment, diagnostics, phylogeny) for the final call.

If the goal is “define SINE families,” neither algorithm is the authority.  
If the goal is “propose and audit structure in SINE space,” HDBSCAN is often worth running — especially beside, not instead of, vsearch.
