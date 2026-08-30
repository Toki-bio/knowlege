---
title: Coalescence as same chromosome copy
type: explanation
---

> Type: explanation

# Coalescence is two DNA lineages becoming one copy, not two people becoming one ancestor

Pedigree thinking and coalescent thinking use the word “ancestor” for different objects.

A pedigree ancestor is a **person**.  
A coalescent ancestor is a **chromosome copy**: a physical DNA molecule that was replicated and passed on.

Until that distinction is sharp, “time to the most recent common ancestor,” heterozygosity, incomplete lineage sorting, and PSMC all feel like word games.

---

## A chromosome copy is a lineage of molecules

Ignore recombination for a few sections. Follow one autosome, say chromosome 5.

A woman living 100,000 years ago had two chromosome-5 molecules:

```text
copy A  ←  the chromosome 5 she got from her mother
copy B  ←  the chromosome 5 she got from her father
```

Those are two different pieces of DNA, even though they sit in one nucleus.

When she has a child, that child inherits **one** of them (still ignoring recombination). Suppose the child gets copy A. A grandchild can inherit that same A lineage, then a great-grandchild, and so on. Each generation the molecule is replicated; the lineage is still “copy A.”

So “the same chromosome copy” does not mean the same physical atoms surviving for millennia. It means an unbroken **copying lineage**.

---

## Walk your two copies backward

You have two chromosome-5 molecules today:

```text
🔴  maternal copy
🔵  paternal copy
```

Each generation, looking backward, each lineage has to pick **exactly one** parental molecule to continue through. The history looks like two threads:

```text
today
🔴                🔵
 │                │
 │                │
 │                │
```

Sooner or later those threads occupy the **same** ancestral molecule:

```text
today
🔴                🔵
 │                │
 └───┐        ┌───┘
     │        │
     ─── 🔴 ───     ← one ancestral copy
```

From that generation farther back there is only one thread left. **That merger is coalescence.** The person who carried that molecule is incidental. What merged is the DNA lineage.

The most recent such molecule is the **most recent common ancestral copy** of the two sequences. Its age is the local TMRCA.

---

## Same person is not coalescence

Suppose both of your threads, walking backward, enter the same woman 200,000 years ago. Have they coalesced?

Not yet, if they still sit on her two different chromosome-5 copies: maternal lineage on copy A, paternal lineage on copy B. She is a shared **pedigree** ancestor. The DNA lineages are still two molecules, and they can still pick up different mutations.

Coalescence happens only when both threads enter **the same copy**.

This is why coalescent models are written as if the population were haploid. A diploid population of *N* people is treated as **2*N* gene copies**. Rosenberg & Nordborg (2002, Fig. 4) draw exactly that: a genealogy of copies, then note that diploids are handled by doubling the haploid pool. The “parent” a lineage picks in that figure is a parental **gene copy**, not a person with two chromosomes.

---

## Why the definition is about copies: mutations

Mutations happen on molecules.

While 🔴 and 🔵 are separate lineages, each can mutate. Those mutations become differences between the two sequences you hold today (heterozygous sites, if you are looking inside one person).

Once the two threads have coalesced, there is only one ancestral molecule left. Mutations on that shared lineage are inherited by **both** descendants. They do not distinguish 🔴 from 🔵.

Looking forward from the common copy: the split creates two lineages; differences between today’s sequences must have arisen on those two branches. Looking backward: after coalescence you cannot generate further differences between the pair, because there is only one molecule to mutate.

That is the whole bridge:

```text
mutations happen on copies
        ↓
differences exist only while copies have separate histories
        ↓
the age of the shared copy (TMRCA) sets how much difference you expect
        ↓
the distribution of those ages across the genome is a record of population history
```

A large population has many parental copies to pick from, so two sampled lineages wait longer to hit the same one. Coalescence is slow; TMRCAs are old; more mutations accumulate; diversity is high. A small population does the opposite. Effective population size *N*e is, in this picture, “how many copies were there to miss each other.”

---

## Recombination: one chromosome, many coalescences

The story above is for one non-recombining piece of DNA. Real autosomes recombine.

A recombination event in an ancestor splices two different parental molecules into one descendant chromosome. Walking backward, that is a **split** of the ancestral thread: the left part of the sequence continues through one parental copy, the right part through another.

So neighboring segments of the same chromosome can have **different** copying histories and different TMRCAs:

```text
          recombination in some ancestor
                    │
     left segment   │   right segment
     coalesces      │   coalesces
     deep           │   recently
     (old TMRCA)    │   (young TMRCA)
```

A diploid genome is therefore a mosaic of local genealogies. Heterozygous stretches tend to sit on segments whose two copies have been separate for a long time (old local TMRCA, more chance to mutate). Long homozygous stretches tend to sit on segments that coalesced recently.

The full object that records every coalescence and every recombination in a sample is an **ancestral recombination graph** (ARG). PSMC is a stripped-down reading of the two-copy slice of that object.

---

## What PSMC is doing with that mosaic

[PSMC](https://github.com/lh3/psmc) (Li & Durbin 2011) takes **one diploid genome**. It does not need a pedigree. It needs a map of which sites are heterozygous.

It treats the genome as a hidden Markov path along the chromosome. The hidden state at each window is the local coalescence time of the two copies. The observations are “this 100 bp bin had a heterozygote” vs “it did not.” Mutation rate turns expected heterozygosity into a guess at TMRCA; recombination rate turns neighboring windows’ TMRCAs into a guess at where the genealogy switched.

Then it uses the **density of those TMRCAs through time** to infer *N*e(*t*): epochs with few coalescences look like large populations; epochs with many coalescences look like small ones.

Mather, Traves & Ho (2020) Fig. 1a is the picture of this: two present-day copies, local trees that are just a single internal branch (two leaves have only one topology), recombination breakpoints, and a different TMRCA in each segment. Those TMRCAs *are* PSMC’s hidden states.

PSMC is weak very recently (too few recombination breakpoints to time young coalescences well) and assumes a single panmictic population. Structure, inbreeding, and bad heterozygote calls all masquerade as size change. The conceptual point still stands: the method is reading a mosaic of **copy coalescences**, not a family tree of people.

---

## One continuous picture

```text
a mutation arises on one chromosome copy
        │
that copy is replicated down some descendants
        │
today you hold two copies (🔴 / 🔵)
        │
walk both backward; each generation each thread
picks one parental molecule
        │
they enter the SAME molecule  =  coalescence
(same person with two still-separate copies is not enough)
        │
mutations after that split  →  differences between 🔴 and 🔵
mutations on the shared stem →  shared by both, not differences
        │
recombination cuts the chromosome into segments
with different coalescence times
        │
today: mosaic of heterozygous / homozygous stretches
        │
PSMC (and ARG methods) infer local TMRCAs from that mosaic,
then Ne(t) from how often coalescences fell in each epoch
```

If a diagram only shows people merging into ancestors, it has skipped the molecule.

---

## Figures that actually show this

These are the pictures that match the argument above, not a leftover reading list.

1. **Rosenberg & Nordborg (2002).** *Genealogical trees, coalescent theory and the analysis of genetic polymorphisms.* *Nature Reviews Genetics* 3: 380–390. [doi:10.1038/nrg795](https://doi.org/10.1038/nrg795)  
   Fig. 1: polymorphism as mutations on a genealogy of copies. Fig. 2: random coalescent trees (same model, different chance outcomes). Fig. 4: sampled copy lineages walking back through a population of copies; diploids = 2*N* copies.

2. **Mather, Traves & Ho (2020).** *A practical introduction to sequentially Markovian coalescent methods.* *Ecology and Evolution* 10: 579–589. [PMC7045566](https://pmc.ncbi.nlm.nih.gov/articles/PMC7045566/)  
   Fig. 1a: the two-copy, along-the-genome view that PSMC uses.

3. **Rasmussen, Hubisz, Gronau & Siepel (2014).** *Genome-wide inference of ancestral recombination graphs.* *PLOS Genetics* 10: e1004342. [doi:10.1371/journal.pgen.1004342](https://doi.org/10.1371/journal.pgen.1004342)  
   Fig. 1: chromosomes at the bottom, lineages up through time, recombination and coalescence marked, different local trees for different segments. More than two copies, but the clearest ARG cartoon.

4. **Li & Durbin (2011).** *Inference of human population history from individual whole-genome sequences.* *Nature* 475: 493–496. [doi:10.1038/nature10231](https://doi.org/10.1038/nature10231)  
   The original PSMC: one diploid genome → local TMRCAs → *N*e through time.

Informal diagrams of the same two-copy story also appear in various “watchmaker / practical pop-gen” notes; they are teaching aids, not the source of the model.

---

## Takeaways

* Coalescence merges **gene-copy lineages**. Meeting the same ancestral person is not enough if the copies are still two molecules.
* Differences between sequences are mutations that happened **while the copies were separate**. The shared stem cannot create those differences.
* Recombination makes a chromosome a mosaic of local TMRCAs. Heterozygosity is a noisy map of that mosaic.
* PSMC reads the mosaic in one diploid genome. It is coalescent inference on copies, not pedigree reconstruction.
* Standard theory already treats diploids as a pool of 2*N* copies. The “people vs DNA” confusion is what that modeling choice is there to avoid.
