---
title: Newmap, HiTE, and NSIT for SINE search
type: explanation
---

> Type: explanation

# Newmap, HiTE, and NSIT: what each tool is actually for

## The starting question

Newmap is an effective tool for finding **unique k-mers / uniquely mappable read lengths** in a genome.

Can it be applied to find **non-perfect matches**, such as novel or divergent SINEs?

**Short answer: no — not as a SINE searcher.**

Newmap answers a different question. Confusing mappability tools with homology search is how this topic usually goes wrong.

---

## What Newmap actually does

[Newmap](https://github.com/hoffmangroup/newmap) is a **genome mappability** package.

For each position in an assembly, it asks roughly:

> How long does a read starting here have to be before that sequence occurs only once in the genome?

From those minimum unique lengths it can produce:

* single-read mappability tracks,
* multi-read mappability tracks.

It is excellent when you need to know which genomic intervals are uniquely placeable at a given read length.

### What it does not do

* It does **not** search for approximate matches.
* It does **not** detect sequences that are similar-but-not-identical to a query.
* It does **not** discover novel SINEs by allowing mismatches / indels relative to a known consensus.

So if the biological goal is:

> Find SINE-like sequences that diverge from known models, or insertions absent from a reference,

Newmap is the wrong class of tool.

---

## Why novel / divergent SINEs need approximate matching

SINE detection is rarely an exact-k-mer problem.

Real copies differ from consensuses by:

* substitutions,
* small indels,
* truncated ends,
* lineage-specific diagnostic changes,
* degradation in old copies.

Exact unique k-mers are useful for **mappability** and for some seeding strategies, but uniqueness is almost the opposite of what you want when hunting a repetitive family: SINEs are interesting precisely because related copies exist in many places.

The right question is usually:

> Which genomic intervals are similar enough to my SINE model / novel insertion to count as members?

not:

> Which k-mers are unique in the genome?

---

## Better tool classes for non-perfect SINE detection

| Goal | Better class of tool | Notes |
| --- | --- | --- |
| Unique exact k-mers / mappability | **Newmap** | Exact uniqueness only |
| Unique k-mers with small edit distance | tools like UniqueKMER / related filters | Still about uniqueness, not family search |
| Approximate mappability with mismatches | **GenMap**-style `(k,e)` mappability | Allows mismatches, still mappability-oriented |
| Sequences absent from a reference assembly | **NSIT**-style novel-sequence detection | Finds novelty relative to a reference, not “all SINEs” |
| De novo full-length TE library from an assembly | **HiTE**, RepeatModeler2, EDTA, etc. | Builds / discovers TE models |
| Search a known SINE consensus across a genome | homology search: **ssearch / BLAST / nhmmer / RepeatMasker `-lib`** | This is library search, not mappability |
| Sensitive profile search for degenerate copies | **nhmmer** / profile HMMs | Strong for old / divergent copies |

The important separation:

```text
mappability tools     →  "Is this position uniquely placeable?"
novel-sequence tools  →  "What is absent from the reference?"
de novo TE tools      →  "What TE models exist in this assembly?"
homology / HMM search →  "Where are copies of this SINE family?"
```

Novel SINE work often needs the last three. Newmap belongs to the first.

---

## NSIT: useful for novelty, not for all SINE copies

[NSIT](https://doi.org/10.1371/journal.pone.0108011) (Novel Sequence Identification Tool) was designed to find sequences present in a de novo assembly but **absent from a reference assembly**.

That can help when the question is:

> Which insertions / sequence blocks are new relative to my reference?

Those novel blocks may include young SINE insertions, but NSIT is not a SINE annotator. It does not by itself:

* classify SINE structure,
* demultiplex subfamilies,
* or recover old degenerate copies already present in the reference.

So NSIT is a **candidate generator for novelty**, not a complete SINE discovery pipeline.

---

## HiTE: a fast TE detector — but not “RepeatMasker with my SINE library”

[HiTE](https://github.com/CSU-KangHu/HiTE) ([Hu et al., 2024](https://doi.org/10.1038/s41467-024-49912-8)) is a strong tool for **de novo full-length TE detection and annotation in genome assemblies**.

Why people prefer it to RepeatMasker-centered workflows:

* RepeatMasker with a large TE library on a mammalian genome can be painfully slow.
* HiTE is designed to discover structurally intact TE models efficiently and can be much faster than older de novo library-building stacks in published comparisons.
* It can then annotate the genome with the library it generated (`--annotate 1`), producing coordinate outputs such as `HiTE.out` / `HiTE.gff`.

### Important correction to a common ChatGPT claim

It is easy to describe HiTE as:

> “Feed in my SINE consensus library and get all approximate genome matches, like a faster RepeatMasker.”

That is **not** the accurate primary use case.

In HiTE’s documented interface:

* `--genome` is required — it works on an assembly;
* `--curated_lib` is a trusted library used to **pre-mask** highly homologous sequence and reduce work;
* `--annotate` annotates using the **TE library HiTE generated**.

So HiTE is best thought of as:

```text
assembly
   ↓
de novo full-length TE discovery
   ↓
TE library (confident_TE.cons.fa, etc.)
   ↓
optional genome annotation
```

not as a general-purpose “query my custom SINE FASTA against genome” replacement for `RepeatMasker -lib`, `nhmmer`, or `ssearch36`.

If you already have a SINE consensus and only want genome-wide approximate copies, a targeted homology / HMM search is usually the direct tool. HiTE is more useful when you want a **fast de novo TE library / full-length TE annotation** from the assembly itself.

---

## A sane hybrid, with roles kept honest

The useful idea from the exchange is still valid if the tools are assigned the right jobs:

```text
1. Novelty / candidate generation
   NSIT (or equivalent): sequences absent from reference
        ↓
2. SINE-like filtering
   length, poly-A / 3′ tail, internal structure, homology to known SINEs
        ↓
3. Collapse and consensus
   vsearch / cd-hit-est → consensus or subfamily consensuses
        ↓
4a. Targeted genome search for that family
   ssearch / BLAST / nhmmer / RepeatMasker -lib with ONLY that library
        ↓
4b. Parallel genome-wide TE discovery (optional)
   HiTE on the assembly for full-length TE models / annotation
        ↓
5. Optional refinement
   build HMM from accepted copies → nhmmer for ultra-degenerate members
```

Why this is often better than “just run RepeatMasker on everything”:

* you control the search library when looking for a specific novel SINE;
* you avoid scanning an enormous generic TE library when that is not the question;
* HiTE can still be used where it shines: building / annotating full-length TE content of the assembly;
* Newmap can still be used later if you care about mappability of the resulting loci — but not as the discovery engine.

---

## Where Newmap still belongs in a SINE project

Newmap is not useless here. It is just downstream / orthogonal.

Useful roles:

* ask whether a candidate “novel SINE insertion” lands in uniquely mappable sequence;
* interpret short-read support around a SINE locus;
* avoid over-interpreting pileups in low-mappability flanks;
* compare mappability of young insertions vs ancient degenerate copies.

Useless roles:

* discovering divergent SINE copies by approximate match;
* replacing RepeatMasker / nhmmer / ssearch;
* clustering SINE subfamilies.

---

## Practical decision guide

**If you want unique mappability tracks**  
→ Newmap (exact) or GenMap-style `(k,e)` tools (with mismatches).

**If you want sequence absent from a reference**  
→ NSIT-style novel-sequence detection, then filter for SINE-like structure.

**If you want a de novo full-length TE library from an assembly quickly**  
→ HiTE (and compare with RepeatModeler2 / EDTA if needed).

**If you want all approximate copies of a specific SINE consensus**  
→ targeted homology / profile search (`ssearch`, BLAST, `nhmmer`, or RepeatMasker with a small custom `-lib`), not Newmap and not “HiTE as a query searcher.”

**If RepeatMasker is too slow**  
→ first ask *why*: huge generic library? whole-genome annotation? custom family search?  
Then replace the slow step with the matching tool, instead of assuming one program solves every TE problem.

---

## Key takeaways

* Newmap finds uniquely mappable exact sequence contexts. It does not find novel / divergent SINEs by approximate matching.
* Novel SINE detection needs homology, novelty detection, and/or de novo TE modelling — not uniqueness alone.
* NSIT is for sequences absent from a reference; HiTE is for fast de novo full-length TE detection/annotation in assemblies.
* Do not describe HiTE as a drop-in “custom SINE library searcher” equivalent to RepeatMasker `-lib`.
* The productive hybrid is: novelty/filter/consensus with NSIT-like methods, targeted search for that family, and HiTE when you need assembly-wide TE discovery — with Newmap reserved for mappability questions.
