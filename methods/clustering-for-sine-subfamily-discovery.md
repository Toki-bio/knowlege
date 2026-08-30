---
title: Clustering for SINE subfamily discovery
type: explanation
---

> Type: explanation

# Clustering for de novo SINE subfamily discovery

## The question this note answers

Could **K-means** or **Gaussian Mixture Models (GMMs)** help with **de novo SINE discovery and classification** — for example after a pipeline like [SINE-de-novo-genome-scan](https://github.com/Toki-bio/SINE-de-novo-genome-scan)?

Short answer:

> They are **not** established methods for assigning SINE subfamilies.  
> They **are** a plausible next step for **candidate discovery**: proposing groups that still need biological validation.

That distinction matters. Treating a cluster label as a subfamily name is how clustering becomes cargo-cult taxonomy. Treating a cluster as a hypothesis is how it becomes useful.

---

## What "subfamily" already means in SINE work

In practice, SINE subfamilies are not defined by an unsupervised algorithm. They are defined by a combination of:

* shared diagnostic mutations relative to a parent consensus,
* recoverable consensus sequences,
* phylogenetic coherence (often monophyly among copies),
* structural hallmarks (A-box / B-box promoters, poly(A) or related 3′ tails, target-site duplications),
* and, when available, placement against curated libraries (Dfam, Repbase / RepeatMasker libraries, SINEbase, manuscript-specific consensuses).

The standard toolkit is therefore homology- and model-driven:

| Approach | Typical role |
| --- | --- |
| BLAST / ssearch / RepeatMasker | detect copies against known models |
| Consensus reconstruction | define the unit being classified |
| Profile HMMs | sensitive recognition of known families |
| Phylogenetics | test whether a proposed group hangs together |
| Curated databases | name and compare against prior annotation |
| Supervised ML (increasingly) | classify into already-labelled classes |

K-means and GMMs do not replace that stack. If they enter the workflow at all, they enter as a **proposal engine** sitting between "I have many SINE-like copies" and "I claim a new subfamily."

---

## Two kinds of clustering that are easy to confuse

SINE pipelines already "cluster" in several places. Those steps are not the same as GMM-based subfamily discovery.

### 1. Sequence-identity clustering (redundancy reduction)

In `SINE-de-novo-genome-scan`, fragments of known SINE consensuses are clustered with `vsearch` at **85% identity** before genome search. That step asks:

> Which query fragments are near-duplicates of each other?

Its job is to thin a search library, not to invent subfamilies in a new genome.

### 2. Genomic hit merging (locus assembly)

Later, genome hits within a window (e.g. 500 bp) are merged into candidate loci. That step asks:

> Which hits probably belong to the same insertion?

Again: locus definition, not subfamily taxonomy.

### 3. Feature-space clustering (the proposal discussed here)

This is different. After you have a set of high-confidence SINE copies from a genome, you represent each copy as a **numeric feature vector** and cluster in that space. That step asks:

> Do these copies fall into natural groups that might correspond to undescribed lineages?

Only this third sense is what K-means / GMM are being considered for.

```
Sequence-identity clustering   →  reduce redundant queries
Genomic hit merging            →  assemble loci
Feature-space clustering       →  propose candidate subfamilies
```

---

## Where clustering belongs in a discovery pipeline

A useful placement is **after detection, before naming**:

```
Genome
   ↓
SINE detector
(e.g. fragment search → locus merge → candidates.fa)
   ↓
Quality filter
(drop truncated / low-confidence copies)
   ↓
Feature extraction
   ↓
Clustering (K-means / GMM / Bayesian mixture)
   ↓
Candidate groups
   ↓
Biological validation
(consensus, diagnostics, phylogeny, structure)
   ↓
Subfamily proposal  — or rejection
```

The detector remains responsible for "is this a SINE-like insertion?"  
Clustering is responsible only for "do these insertions look like they form separable groups?"

That matches the spirit of existing de novo tribe / young-lineage clustering in other SINE workflows: unsupervised grouping is a starting point, not the final taxonomy.

---

## Why the method is hypothesis generation, not classification

A cluster is cheap. A subfamily is expensive.

For a cluster to graduate into a candidate subfamily, it should survive checks such as:

* **Shared diagnostics** — do members share substitutions or indels that distinguish them from neighbouring groups?
* **Consensus recoverability** — does a stable consensus emerge, or is the cluster a grab-bag?
* **Phylogenetic coherence** — do the copies form a clade, or at least a coherent cloud, relative to related SINEs?
* **Promoter conservation** — are A-box / B-box (or family-appropriate) motifs coherent within the group?
* **Age and abundance** — is the group young and expanding, old and eroded, or an artifact of fragmentation?
* **Genomic distribution** — does it behave like a real mobile-element lineage rather than a local assembly quirk?
* **Independence from known models** — is it truly novel, or just a degraded / chimeric version of something already named?

Only after those checks does a cluster become biologically interesting.

---

## K-means versus GMM

### K-means

K-means partitions points into *K* groups by minimizing within-cluster squared distance to centroids.

**Strengths**

* simple and fast
* scales to large copy counts
* good first look at structure in a feature space

**Weaknesses for SINE data**

* hard assignments: every copy is forced into exactly one cluster
* assumes roughly compact, spherical clusters
* struggles with transitional / mosaic sequences
* requires choosing *K* up front

K-means is a useful baseline. It is a poor final authority on evolutionary boundaries.

### Gaussian Mixture Models

A GMM models the data as a weighted sum of Gaussians. Each copy gets a **probability** of belonging to each component.

**Strengths**

* soft assignments
* overlapping groups are allowed
* ambiguous copies remain visible instead of being forcibly labelled

**Weaknesses**

* still assumes Gaussian-shaped components in feature space
* sensitive to how features are scaled and chosen
* still needs a principled way to choose the number of components (unless made Bayesian)

For SINE evolution, soft assignment is the biologically attractive part.

---

## Why soft membership matters

SINE subfamilies often arise by progressive accumulation of changes:

```
Subfamily A
      \
       \
     transitional copies
            \
             \
           Subfamily B
```

Hard clustering forces a false choice:

```
A
or
B
```

A GMM can instead say:

```
P(A) = 0.63
P(B) = 0.37
```

Those intermediate copies are not noise to be cleaned away. They are often the most interesting objects in the dataset: possible nascent lineages, incompletely sorted variants, or chimeras.

In other words, GMM is useful not only for the clusters it finds, but for the **uncertainty** it preserves.

---

## What to cluster on

Raw k-mer frequencies are an easy default and usually a weak one. They mix true diagnostic signal with composition noise, truncated ends, and homology to unrelated AT-rich sequence.

More useful is a compact feature vector built from biology-facing quantities, for example:

| Feature class | Examples |
| --- | --- |
| Divergence | distance to nearest known consensus; private substitutions |
| Diagnostics | presence/absence of candidate defining mutations |
| Indels | characteristic insertions/deletions relative to a parent model |
| Structure | A-box / B-box integrity; tail length and purity |
| Flanks | TSD length/identity; insertion motif |
| Age proxies | substitution level; CpG decay where relevant |
| Embeddings | learned sequence vectors, if used carefully and validated |

The representation decides the biology you are allowed to discover. Cluster identity-space and you rediscover percent-identity tribes. Cluster diagnostic-mutation space and you are closer to classical subfamily logic.

---

## Small clusters and outliers are the point

In a real genome scan, the largest clusters are usually boring in the best sense: they recover abundant, already-known or easily consensus-able groups.

The scientifically valuable residue is often:

* small but stable clusters,
* copies with middling probability under every major component,
* outliers poorly explained by the whole mixture.

These may be:

* previously undescribed subfamilies,
* lineage-specific variants,
* chimeric or rearranged SINEs,
* or simply assembly / detection artifacts.

So the analysis should be designed to **surface** rare structure, not only to paint the dominant clouds.

---

## Why Bayesian mixtures are attractive

Ordinary GMM still asks you to choose *K*. For de novo SINE work, that is awkward because the number of real lineages is exactly what you do not know.

A **Bayesian Gaussian mixture** (or a Dirichlet-process mixture) is attractive because:

* the number of components can be inferred rather than fixed by hand,
* unused components are driven toward vanishing weight,
* uncertainty about *K* becomes part of the model instead of a silent tuning knob.

That does not remove the need for biological validation. It only reduces one arbitrary choice in the proposal stage.

Conceptual target:

```
High-confidence SINE copies
        ↓
Feature matrix
        ↓
Bayesian / Dirichlet-process GMM
        ↓
Candidate components + membership probabilities
        ↓
Consensus + diagnostics + phylogeny
        ↓
Keep / split / reject each component
```

---

## Concrete placement after SINE-de-novo-genome-scan

The public pipeline already ends near the right hand-off point. After `sine_scan.sh` you have:

* `candidate_loci.bed`
* `candidates.fa`

That is the natural insertion point for feature-space clustering:

```
sine_scan.sh
   → candidates.fa
        ↓
filter to full-length / high-quality copies
        ↓
align to nearest seed consensus (or all-vs-all among candidates)
        ↓
build feature table
        ↓
Bayesian GMM
        ↓
per-cluster consensuses
        ↓
SubFam / SINEderella-style assignment logic
   or phylogenetic tests against known libraries
```

Important boundary:

* **Do not** replace the detector with clustering.
* **Do not** treat GMM labels as final subfamily names.
* **Do** use clustering to prioritize which candidate sets deserve consensus reconstruction and manual/phylogenetic review.

This also clarifies what the existing `vsearch --id 0.85` step is *not*: it is library compression before search, not subfamily discovery in the target genome.

---

## Practical failure modes

If this approach is tried carelessly, the usual failure modes are predictable:

1. **Clustering truncated junk** — 5′ or 3′ fragments form "subfamilies" that are really incompleteness classes.
2. **Clustering composition** — AT-rich tails and low-complexity sequence dominate the metric.
3. **Re-discovering detection bias** — copies found by the same seed fragment cluster together because of search ancestry, not biology.
4. **Over-splitting by age** — old degraded copies of one lineage shatter into many clusters.
5. **Under-splitting by homology** — distinct young subfamilies remain too close in feature space and collapse into one component.
6. **Naming too early** — exporting cluster IDs as if they were curated subfamily names.

Mitigations follow directly: filter before clustering, choose features that encode diagnostics rather than raw composition, validate every kept component with consensus and phylogeny, and keep cluster IDs internal until validation succeeds.

---

## What this is not claiming

To keep the claim honest:

* K-means / GMM are **not** presented here as standard SINE subfamily classifiers in the literature.
* This note does **not** assert that published SINE taxonomies were built primarily by mixture models.
* Integration with `SINE-de-novo-genome-scan` is architectural: the natural hand-off is after `candidates.fa`, not a claim that the repository already implements GMM subfamily calling.

The valuable idea from the discussion is narrower and stronger:

> After reliable SINE detection, unsupervised mixture modelling in a biologically chosen feature space is a reasonable way to **propose** candidate subfamilies — especially when assignments are probabilistic and every proposal is forced through consensus and phylogenetic validation.

---

## Key takeaways

* Use homology, consensuses, HMMs, phylogenetics, and curated libraries to define SINE subfamilies.
* Use K-means / GMM only as a **hypothesis generator** after high-confidence copies exist.
* Distinguish redundancy clustering and locus merging from feature-space subfamily proposal.
* Prefer soft, preferably Bayesian, mixtures over hard K-means when evolutionary boundaries are gradual.
* Feature choice is the biology; the clustering algorithm is only the grouping machinery.
* Small clusters, ambiguous memberships, and outliers are often more informative than the giant known groups.
* No cluster becomes a subfamily until diagnostics, consensus, and phylogeny agree.
