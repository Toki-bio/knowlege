---
title: Coalescence as same chromosome copy
type: explanation
---

> Type: explanation

# Coalescence is about DNA molecules, not people

This is exactly where almost everyone gets stuck: we naturally think about **people**, while coalescent theory thinks about **DNA molecules**.

## One sentence that often makes it click

Do **not** think of coalescence as:

> Two ancestors becoming one.

Think of it as:

> Two photocopies eventually tracing back to the same original document.

Today you have two photocopies (maternal and paternal chromosome copies). Trace their copying history backward and eventually both were copied from the same original page. From that point farther back, there is only one page to follow.

That “same original page” is what geneticists mean by the **same ancestral chromosome copy**: a lineage of DNA, not merely the same ancestral person.

---

## Forget people for a moment — follow one DNA molecule

Suppose 100,000 years ago there was a woman.

She had two copies of chromosome 5:

```text
copy A  ← chromosome 5 from her mother
copy B  ← chromosome 5 from her father
```

These are two different physical DNA molecules.

She has a child. Ignoring recombination for the moment, the child inherits **only one** of these copies.

Suppose the child gets copy A. Then a grandchild may inherit that same A lineage. Generation after generation, this can be the **same chromosome copy lineage**, simply being replicated.

---

## Now follow your two chromosomes backward

Today you have:

```text
🔴 maternal chromosome 5
🔵 paternal chromosome 5
```

Trace them backward. Every generation, each lineage chooses exactly one parental chromosome copy to follow:

```text
Today

🔴          🔵
│           │
│           │
│           │
```

Eventually something like this happens:

```text
Today

🔴          🔵
│           │
│           │
└──────┐ ┌──┘
       │ │
       🔴
```

Suddenly both lineages are following the **very same ancestral chromosome copy**.

Not just the same person.  
Not just the same ancestor.  
The **same physical lineage of DNA**.

From that point backward there is only one lineage left.

**That is coalescence.**

---

## Why “same person” is not enough

Suppose 200,000 years ago both of your lineages reach the same woman.

Did they coalesce?

**Not necessarily.**

That woman has two chromosome-5 copies.

* your maternal lineage might enter one copy,
* your paternal lineage might enter the other.

They are still different DNA molecules.

Only when both lineages enter the **same chromosome copy** has coalescence occurred.

This is why coalescent theory tracks **gene copies / chromosome copies**, not individuals.

---

## Why this matters for mutations

Mutations occur on chromosome copies.

**Before coalescence:**

* one lineage accumulates mutations,
* the other lineage accumulates mutations,
* independently.

**After coalescence (looking backward):**

* there is only one ancestral lineage left,
* no further differences can arise between the two present-day copies *before* that point, because there is only one DNA lineage to mutate.

Looking forward from the common ancestral copy: every mutation on that shared lineage is inherited by both descendant copies, so it does not create a difference between them. Differences between today’s 🔴 and 🔵 must have arisen on the separate branches after they split.

That is the bridge from mutations → coalescence → genetic diversity.

---

## Recombination makes many local histories

The story above ignored recombination.

With recombination, different segments of a chromosome can have **different genealogies**:

```text
left segment coalesces long ago
right segment coalesces more recently
```

So a diploid genome is a mosaic of local TMRCAs (times to most recent common ancestral copy). Methods like **PSMC** exploit that mosaic: local heterozygosity patterns reflect local coalescence times, which in turn inform historical effective population size.

Coalescence is still “same ancestral DNA copy” — just applied locus by locus along the chromosome.

---

## Good graphical sources

These figures are especially useful:

### 1. Best for PSMC specifically

**Mather et al. (2020)** — *A practical introduction to sequentially Markovian coalescent methods*  
[PMC article](https://pmc.ncbi.nlm.nih.gov/articles/PMC7045566/)

Look at **Figure 1a**.

It shows what PSMC-style models use:

* two present-day chromosome copies,
* local genealogies,
* recombination breakpoints,
* different TMRCA for neighboring regions.

Almost a graphical explanation of the PSMC hidden states.

### 2. Best for “same ancestral chromosome copy”

**Nature Reviews Genetics** review on genealogical trees / coalescent theory / polymorphisms  
(search: *Genealogical trees, coalescent theory and the analysis of genetic polymorphisms*)

Especially useful:

* Figure 2 — random genealogical trees
* Figure 4 — basic principle of coalescence

These make it easier to see lineages merging backward in time as **gene copies**, not people.

### 3. Best for recombination changing genealogy

**Genome-Wide Inference of Ancestral Recombination Graphs** (*PLOS Genetics*)

**Figure 1** is outstanding: chromosomes at the bottom, lineages traced backward, recombination events, coalescence events, and different trees for different segments.

Even with more than two chromosomes, it is one of the clearest pictures of why neighboring loci can have different ancestral histories.

### 4. Surprisingly good intuitive overview

**The Watchmaker’s Guide to Population Genetics** — PSMC overview pages

Not primary literature, but often has very readable diagrams of two chromosome copies, MRCA, heterozygous sites, recombination, and PSMC inference.

---

## A useful custom schematic to keep in mind

One continuous mental diagram that connects the pieces:

```text
1. A mutation arises on one chromosome copy
        ↓
2. That copy is inherited through generations
        ↓
3. Today: two chromosome copies (🔴 / 🔵)
        ↓
4. Trace both backward until they occupy
   the same ancestral chromosome copy
   = coalescence
        ↓
5. Recombination splits the chromosome into
   segments with different coalescence times
        ↓
6. Today you observe a mosaic of
   heterozygous / homozygous stretches
        ↓
7. PSMC-like methods infer historical Ne
   from that mosaic of local TMRCAs
```

If a figure does only “people merging into ancestors,” it has skipped the essential step.

---

## Key takeaways

* Coalescence merges **gene / chromosome copy lineages**, not merely people.
* Meeting the same ancestral person is not enough if the lineages still sit on different chromosome copies.
* Mutations create differences only while lineages are separate; after coalescence (backward), there is one ancestral DNA lineage.
* Recombination creates a mosaic of local coalescence times along a chromosome.
* PSMC uses that mosaic; good intuition figures are Mather et al. Fig. 1a, classic coalescent review figures, and ARG Figure 1.
