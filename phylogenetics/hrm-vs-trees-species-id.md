---
title: HRM vs trees for species ID
type: explanation
---

> Type: explanation

# Why an HRM paper does not build a phylogenetic tree

Paper context: Distinguishing fanged frogs (*Limnonectes*) from Thailand with high-resolution melting analysis — a diagnostic assay for the *Limnonectes kuhlii* complex using mitochondrial **16S**.

## Short answer

They do not build a tree because the paper is **not trying to infer evolutionary relationships**.

It is trying to build a **cheap, fast diagnostic test** for known species. Melting curves are good for that. Trees are for different questions.

---

## What they actually do

1. PCR-amplify a short fragment of mitochondrial **16S rRNA**.
2. Gradually heat the amplicon after PCR.
3. Measure fluorescence as the duplex melts.
4. Read a **melting curve** and melting temperature (**Tm**).

Different sequences tend to have different GC content and base-stacking, so they produce different melting profiles.

Instead of sequencing, the assay classifies samples by curve shape / Tm.

Conceptual workflow:

```text
DNA extraction
   ↓
PCR (with saturating dye)
   ↓
HRM melt curve
   ↓
species ID in minutes
```

Versus a sequencing / tree workflow:

```text
DNA extraction
   ↓
PCR
   ↓
sequencing
   ↓
alignment
   ↓
phylogenetic inference
```

The second is slower and more expensive. That is exactly what HRM is designed to avoid.

---

## What HRM is conceptually

Think of HRM as **PCR fingerprinting**, not phylogenetics.

```text
sequence differences → melting curve shape → species label
```

not

```text
sequence → alignment → tree → evolutionary inference
```

Another key point: they are targeting **already known species**. HRM behaves like a barcode classifier, not a discovery method.

It can answer:

> Is this *L. bannaensis* or *L. taylori*?

It cannot answer:

> What is the evolutionary relationship among these frogs?  
> Are there undescribed cryptic lineages here?

---

## Separate three questions that get mixed

| Question | What it asks | HRM? | Tree / sequence? |
| --- | --- | --- | --- |
| **Identification** | Which known label does this sample match? | Yes, if calibrated | Yes, but overkill for routine ID |
| **Species delimitation** | Where are the species boundaries? | No | Much better |
| **Phylogeny** | How are lineages related? | No | Required |

The frog HRM paper targets only **#1**.

Taxonomy / previous sequencing work defines the species.  
The HRM paper builds a quick test for them.

```text
taxonomy paper  → defines species
HRM paper       → builds a rapid diagnostic shortcut
```

---

## If they did a tree, would delimitation be better?

**Yes — if the goal were species delimitation**, a sequence-based tree (or better, multilocus / genomic analysis) would generally be much more informative than HRM.

With sequences and a tree you can examine:

* monophyly,
* genetic divergence,
* coalescent structure,
* support for candidate lineages.

Typical delimitation logic:

```text
samples → alignment → tree / coalescent model
```

Then ask whether individuals form distinct, supported clusters relative to within-group variation. Methods in that space include GMYC, PTP, and related Bayesian / coalescent approaches.

But trees are not automatically perfect either — especially from a single mitochondrial marker such as 16S:

* introgression,
* incomplete lineage sorting,
* mtDNA ≠ species tree.

Modern taxonomy often needs multilocus or genomic data. Still, even a single-gene tree retains far more biological information than a melt curve.

---

## Is melt “perfect enough”?

**No** — not for anything beyond calibrated identification of known taxa.

Your skepticism is the phylogeneticist’s correct instinct.

### What HRM measures

Only physical properties of the PCR product:

* GC content,
* fragment length,
* nearest-neighbor stacking interactions.

Different substitutions can produce nearly identical melting behavior. Different combinations of mutations can cancel out and yield the same curve.

### What a tree / sequence analysis uses

The actual sequence:

* each substitution,
* shared derived characters,
* evolutionary models,
* branch support / coalescent structure.

Information collapse looks like this:

```text
DNA sequence (~200 bp)
        ↓
   alignment
        ↓
      tree
```

versus

```text
DNA sequence
        ↓
   PCR product
        ↓
  one melting curve
        ↓
 1–2 summary features (Tm, curve shape)
```

Most of the sequence information is thrown away.

### When HRM works well

* species already defined,
* marker chosen so known taxa differ enough in melt behavior,
* lab calibrated with authenticated reference samples.

Then it is a diagnostic assay — closer to a medical PCR test than to systematics.

### Why journals still publish these papers

Because the intended users often want applied monitoring:

* wildlife enforcement,
* field labs,
* food authentication,
* disease / pest diagnostics,
* high-throughput biodiversity surveys with known targets.

They want:

```text
sample → PCR → curve → species
```

not

```text
sample → sequence → alignment → tree
```

---

## Honest hierarchy for species work

Rough information content:

```text
genomic phylogeny
   > multilocus phylogeny
      > single-gene tree
         > barcode distances
            > HRM melt curves
```

HRM is near the bottom: useful as a **fast ID assay**, weak as evolutionary evidence.

It cannot reveal:

* cryptic species,
* many misidentifications,
* hybrids,
* ancestry / relationships.

A tree (or better, denser sequence data) can.

This is closely related to the old criticism of early DNA barcoding hype in systematics: classification without evolutionary analysis.

---

## Key takeaways

* The Limnonectes HRM paper skips trees because it is building a rapid identification assay, not studying evolution.
* HRM collapses sequence differences into melt profiles; that is enough for some calibrated ID tasks and insufficient for delimitation or phylogeny.
* Identification ≠ species delimitation ≠ phylogeny — choose the method for the question.
* A tree would be better if the goal were delimitation or discovering cryptic structure; even then, single-mtDNA trees have known limits.
* Melt is not “perfect enough” for evolutionary biology; it is only “good enough” when the species are already known and the assay is carefully validated.
