---
title: Dimensionality reduction for SINE separation
type: explanation
---

> Type: explanation

# Linear vs non-linear dimensionality reduction for SINE separation

## One-sentence summary

**Linear** methods flatten a flat sheet.  
**Non-linear** methods unwrap a crumpled sheet.

For SINE family work:

> Start with **PCA**. Use **UMAP** (or t-SNE) for visualization.  
> Do **not** use t-SNE alone to assign families.

Dimensionality reduction compresses high-dimensional features while trying to keep important structure. It is not a classifier.

---

## Linear dimensionality reduction

Linear methods assume the data roughly lie on a **flat surface** (a line or plane) inside high-dimensional space.

They build new features as **straight-line combinations** of the originals:

```text
new_feature = 2x + 3y − z
```

Intuition: if the cloud is already fairly flat in 3D, you can project it onto a 2D plane without much distortion.

Common methods:

* **PCA** — Principal Component Analysis
* **LDA** — Linear Discriminant Analysis (supervised; needs labels)

**Pros**

* fast and simple
* relatively interpretable
* stable / reproducible
* excellent preprocessing before clustering or modeling

**Limitation**

Fails when the important structure is strongly curved.

---

## Non-linear dimensionality reduction

Non-linear methods assume the data lie on a **curved surface** (a manifold).

They try to preserve neighborhoods or local distances in a flexible way, then unfold the structure into fewer dimensions.

Common methods:

* **t-SNE** — t-Distributed Stochastic Neighbor Embedding
* **UMAP** — Uniform Manifold Approximation and Projection
* **Isomap**

**Pros**

* handles curved / complex patterns
* excellent for visualization
* can reveal local clusters PCA misses

**Limitations**

* slower
* harder to interpret
* sometimes unstable across runs / parameters
* axes and distances are not “real biology units”

---

## Swiss-roll intuition

Imagine points on a rolled sheet:

```text
True structure: a long flat sheet rolled into a spiral
```

### What PCA sees

PCA only understands straight directions. It finds the best flat projection and can squash the spiral into a blob:

```text
Before (true):          After PCA:

   spiral                 squashed cloud
     🌀                      ●●●●●
    🌀🌀         →          ●●●●●●●
   🌀🌀🌀                   ●●●●●●●●
```

Points that were far apart along the unrolled sheet may end up close together. Global curved structure is lost or distorted.

### What UMAP / t-SNE try to do

Non-linear methods attempt to recognize neighborhood structure and unwrap the sheet:

```text
Before (rolled):          After nonlinear:

   🌀🌀🌀                   ●────────●
  🌀🌀🌀🌀       →           ●────────●
 🌀🌀🌀🌀🌀                  ●────────●
```

Result: local neighborhoods are often better preserved, and cluster-like structure becomes easier to see.

Analogy:

* **Linear** ≈ photographing a crumpled map as-is  
* **Non-linear** ≈ unfolding the map before photographing it

---

## When to use which

| Situation | Prefer |
| --- | --- |
| Data roughly flat; want speed / interpretability / preprocessing | **PCA** |
| Complex curved structure; want visualization | **UMAP** (default) or t-SNE |
| Need stable coordinates for downstream modeling | **PCA** first |
| Making a pretty cluster plot | UMAP / t-SNE |
| Defining biological family labels | **neither** — use alignment / similarity / diagnostics |

---

## For SINE family separation specifically

SINE-derived feature spaces (k-mers, alignment scores, divergence metrics, etc.) usually contain:

* some **linear** signal — major divergence, old vs young, large family splits;
* plus **non-linear** structure — subfamilies, drift, messy gradients.

### What each method does for you

#### PCA

* captures major axes of variation;
* often already separates major families;
* stable and reproducible;
* good for sanity checks, initial inspection, and as input to clustering methods such as [HDBSCAN](https://toki-bio.github.io/knowlege/methods/hdbscan-vs-vsearch-for-sine-grouping/).

#### t-SNE

* emphasizes local clusters;
* can make nice-looking islands;
* distances and cluster sizes/spacings are misleading;
* run-to-run variability can be high.

**Do not decide biology from t-SNE alone.**

#### UMAP

* usually a better default visualization than t-SNE;
* preserves local structure and more global structure than t-SNE;
* good for spotting fine substructure PCA misses;
* still not a family assigner.

---

## Critical point for the pipeline

Neither PCA, t-SNE, nor UMAP is a classifier.

They:

* do **not** assign families,
* do **not** provide biological thresholds,
* should **not** be treated as deterministic labels.

For actual SINE separation / assignment, keep using the things that mean something biologically:

* alignment scores / bitscores,
* consensus similarity,
* diagnostic mutations,
* explicit thresholding / voting logic (e.g. SINEderella-style Step 2),
* and human-revisable group boundaries in an [MSA-centered catalog](https://toki-bio.github.io/knowlege/sines/msa-as-sine-consensus-catalog/).

Use dimensionality reduction to:

* validate clusters,
* detect outliers,
* spot mis-assignments,
* decide whether a denser clustering step is warranted.

---

## Practical recommendation

```text
1. Compute features
   (k-mers and/or alignment-derived scores)

2. Run PCA
   → do major families already separate?
   → use as stable coordinates / preprocessing

3. Run UMAP
   → inspect fine structure
   → look for subfamilies PCA flattened

4. Overlay known labels / bitscore assignments
   → mismatches are the interesting cases

5. If needed, cluster in PCA/UMAP space
   (e.g. HDBSCAN)
   → then validate back on sequences / MSA
```

One-line answer for SINE work:

```text
PCA  → backbone, reliable
UMAP → best visualization default
t-SNE → optional, mostly cosmetic
```

---

## Key takeaways

* Linear methods project onto flat subspaces; non-linear methods try to unfold curved manifolds.
* The Swiss-roll cartoon is the right intuition: PCA can squash a spiral; UMAP/t-SNE try to unroll it.
* For SINE families, start with PCA; visualize with UMAP; treat t-SNE as optional cosmetics.
* Dimensionality reduction is for seeing and preprocessing, not for declaring taxonomy.
* Family assignment still belongs to alignment, similarity, diagnostics, and curated thresholds — with clustering / embeddings as audit layers.
