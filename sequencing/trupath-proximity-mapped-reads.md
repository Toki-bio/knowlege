---
title: TruPath proximity-mapped reads
type: explanation
---

> Type: explanation

# TruPath vs mate-pair sequencing: what proximity really buys you

## The question

How does Illumina **TruPath** (TruPath Genome / proximity-mapped reads) generate long-range genomic information from ordinary short-read sequencing, and what does that actually solve for **repetitive DNA**?

The short answer:

> TruPath does not turn short reads into continuous long reads.  
> It preserves evidence that many short reads came from the **same original long DNA molecule**, using **spatial proximity on the flow cell**, and then uses that linkage computationally.

That is powerful for phasing, structural variants, and hard-to-map regions. It is **not** a substitute for spanning an arbitrarily long homogeneous repeat with one continuous sequence.

---

## Classical mate-pair sequencing

In classical mate-pair (or jump) libraries:

1. DNA is fragmented into relatively long inserts, often several kilobases.
2. Fragments are circularized so the two distant ends are brought together.
3. The junction region is sequenced, yielding essentially **two short reads separated by a long insert**.

Conceptually:

```
read -------------------------- read
             ~several kb
```

Mate pairs were historically useful for scaffolding and structural-variant detection. Library preparation was complicated, chimeric junctions were a recurring nuisance, and the approach has largely been superseded by long-read and other long-range methods.

The key property to remember is simple: mate-pair gives you **two linked endpoints** and an approximate insert size. It does not densely sample the intervening DNA.

---

## What TruPath is doing instead

TruPath Genome is Illumina's assay built on **proximity-mapped read technology** (earlier public discussion sometimes called related ideas "constellation" reads). Public descriptions converge on this workflow:

1. Intact DNA — preferably high-molecular-weight, though standard molecular-weight DNA is also used — is introduced onto a patterned flow cell.
2. On-flow-cell **tagmentation** cuts the DNA and attaches adapters while the molecule is spatially captured.
3. Fragments from one original long molecule generate ordinary Illumina clusters.
4. Clusters derived from the same template tend to occur in **neighbouring nanowells**.
5. Downstream software (notably DRAGEN with proximity mode) uses that spatial information to infer which reads are likely linked through the same original molecule.

Conceptually:

```
Mate pair:

read -------------------------- read
             ~5 kb


TruPath / proximity-mapped reads:

read  read  read  read  read  read  read  read
 \     \     |     |     |     /     /     /
           same original long molecule
```

So instead of two ends of one jump, you get **many short Illumina reads distributed along one original molecule**, plus evidence that they belong together.

Illumina's public materials describe reconstruction of long-distance connections on the order of **>200 kb**, with the realized template-length distribution depending strongly on input DNA quality.

---

## Why tagmentation timing matters

Tagmentation simultaneously:

* cuts DNA into sequencing-compatible pieces, and
* attaches Illumina adapter sequences.

In a conventional library, tagmentation or fragmentation happens in a tube. Fragments from many molecules are mixed, and any long-range relationship must be recovered later by barcodes, contiguity-preserving methods, or pairing design.

In TruPath, the critical point is that tagmentation happens **on or at the flow cell after the long molecule has been spatially localized**. Fragments from one molecule therefore remain associated with a local region of the flow cell.

This is why the method should **not** be casually described as "they just barcode every fragment."

---

## How does software know which reads came from the same molecule?

This is the part that should be stated carefully.

Public Illumina / DRAGEN documentation describes a **probabilistic proximity linking model**: reads that are spatially close on the flow cell are assigned linkage probabilities based on spatial and genomic distance. Those probabilities are then reused for mapping, phasing, structural-variant calling, STR recovery, and related analyses.

The safest description is therefore:

> TruPath encodes long-range information as physical/spatial proximity of clusters from the same captured molecule. Computational analysis converts that proximity into probabilistic read-to-read linkage.

What it is **not**, at least in public descriptions:

* not the same as **10x / linked-read barcoding**, where an explicit molecular barcode is written into the sequence of every fragment;
* not a guarantee that every fragment from a molecule is recovered;
* not a continuous long read.

Exact proprietary hardware and inference details should not be over-specified beyond that.

---

## Coverage along the original molecule is patchy

A long molecule can, in principle, be tagmented throughout its length:

```
Original molecule:

|----|----|----|----|----|----|----|----|----|
```

Sequencing is still stochastic. Some fragments never make usable clusters, some are filtered, and local depth varies:

```
██  ███ █    ██ ████   █  ███    ██
```

So a captured ~150 kb molecule is **not** necessarily sequenced continuously from end to end. You obtain many short reads distributed along it, plus proximity evidence connecting those reads.

That is enough to create long-range **context**. It is not the same thing as a HiFi or ONT read that reports one continuous path.

---

## What this does — and does not — do for repeats

This is the crucial limitation, especially for satellite DNA, long tandem arrays, and other low-complexity sequence.

### Case A — repeat flanked by unique sequence

```
unique A ---- REPEAT ---- unique B
```

Reads from inside the repeat may individually map to many places. TruPath does **not** make those repeat bases unique.

But if the same original molecule also carries unique sequence on both sides, ambiguous repeat-derived reads can be associated with the molecule whose unique reads map near A and B.

Proximity therefore supplies **context**, not new sequence information:

```
unique A
   |
   |---- repeat ----|
                    |
                 unique B
```

This is why proximity helps in "dark" or paralogous regions when there is enough informative sequence on the same molecule to anchor placement, phasing, or structural interpretation.

### Case B — molecule wholly inside a huge homogeneous array

```
unique
 |
 |========== multi-megabase satellite array ==========|
```

If the captured molecule is only ~100–200 kb and lies entirely inside the array:

```
       <--------- captured molecule --------->
```

then all of its reads may be repetitive. There is no unique flank telling the software where inside the array the molecule belongs.

Therefore:

**TruPath cannot resolve an arbitrarily long internal repeat merely because the DNA molecule is long.**

It is especially useful when a repetitive or difficult region can be connected to informative sequence elsewhere on the same molecule. It is much less useful when the entire observed molecule is homogeneous repeat.

---

## Comparison with other long-range strategies

A useful mental model:

| Method | What you really get |
| --- | --- |
| **PacBio HiFi** | one long, continuous, highly accurate read |
| **ONT** | one very long read that can span enormous repeats, with higher raw-error characteristics |
| **Mate-pair** | two short reads connected by a known long insert |
| **10x-style linked reads** | many short reads sharing an explicit molecular barcode |
| **TruPath proximity-mapped reads** | many short Illumina reads from one original molecule, linked by flow-cell proximity and probabilistic inference |

So TruPath sits in a distinctive niche:

* keeps Illumina short-read accuracy and SBS economics,
* adds long-range linkage without requiring a continuous long read,
* improves mapping, phasing, SV detection, STR interpretation, and some paralog problems,
* still cannot replace native long-read spanning of giant homogeneous repeats.

Public TruPath / DRAGEN materials emphasize exactly these gains: better placement in difficult regions, megabase-scale phasing with high-molecular-weight input, structural-variant improvements, STR recovery via proximity to flanks, and paralog-aware genotyping — all while remaining short-read data plus linkage.

---

## Practical takeaways for repetitive DNA work

If the biological question is "which haplotype is this variant on?" or "how are these hard-to-map reads connected across tens or hundreds of kilobases?", proximity-mapped reads are highly relevant.

If the biological question is "what is the exact internal sequence of a multi-megabase satellite array?", proximity alone is the wrong tool. You still need a technology that can **span** the array, or independent orthogonal evidence.

In slogan form:

* **Mate-pair:** two ends, one jump.
* **TruPath:** many short reads, one molecule, linked by where they landed on the flow cell.
* **Long reads:** one continuous path through the sequence itself.

Proximity gives you neighbourhood. Continuity gives you the path. For long homogeneous repeats, neighbourhood without an anchor is not enough.

---

## Key takeaways

* TruPath generates long-range information by capturing long DNA, tagmenting it on the flow cell, and using spatial cluster proximity plus computation to link reads from the same molecule.
* It is not classical mate-pair sequencing, and it is not the same as barcode-based linked reads.
* Tagmentation after spatial localization is the enabling step.
* Coverage along a molecule is typically patchy; you do not automatically get continuous end-to-end sequence of the original molecule.
* For repeats, proximity helps when unique or otherwise informative sequence on the same molecule can anchor ambiguous reads.
* If a captured molecule lies entirely inside a long homogeneous repeat, TruPath has no unique place to put it.
* Prefer TruPath-like proximity data for phasing, SV, and difficult mapping; prefer native long reads when the goal is continuous resolution of giant repetitive arrays.
