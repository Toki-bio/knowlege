---
title: mm2-plus vs minimap2
type: explanation
---

> Type: explanation

# mm2-plus vs minimap2: faster parallelism, same aligner

Paper: [Accelerating minimap2 for whole-genome alignment](https://doi.org/10.1093/bioinformatics/btag083) (*Bioinformatics*, 2026)  
Code: [at-cg/mm2-plus](https://github.com/at-cg/mm2-plus)

## Short answer

**No, you probably should not “switch.”**

`mm2-plus` is an optimized implementation of **minimap2**, not a new aligner class. Scientifically, nothing important changes. The rational stance is:

> Keep minimap2 as the default.  
> Use mm2-plus when you want free speed on chromosome-scale assembly-to-assembly jobs.

Don’t fix what isn’t broken. Benchmark once if curious.

---

## What problem the paper actually solves

Modern assemblies often have **few but very long** sequences (chromosome-scale / near-T2T).

Classic minimap2 parallelization often works like:

```text
one thread → one query sequence
```

With only ~20 chromosomes and 48 threads, many CPUs sit idle. Whole-genome alignment is no longer an embarrassingly parallel “many contigs” problem.

So the paper redesigns parallelism and related bottlenecks for that regime.

---

## Real technical novelties

These are engineering / parallel-algorithm changes, not a new biological mapping model.

### 1. Fine-grained parallel chaining

Instead of one thread owning an entire query sequence, mm2-plus parallelizes chaining **inside** a sequence by partitioning anchors (exact matches) so multiple threads can chain one long chromosome at once.

Public description: anchors are partitioned by the **reference sequence they fall on**, then chained independently and aggregated. A chain cannot span more than one reference sequence, so the partitions are independent.

Why this matters:

```text
classic:  24 threads, 10 contigs  → many idle cores
mm2-plus: many threads can work on one long query
```

### 2. Faster primary / secondary chain classification

After chaining, minimap2 spends non-trivial time deciding which chains are primary vs secondary / overlapping. mm2-plus replaces this with a faster **interval-tree-based** classification step.

### 3. Additional engineering speedups

The full set of optimizations also includes:

* parallel sorting in seeding,
* SIMD-accelerated extension / alignment (building on prior mm2-fast-style work),
* SIMD chaining components in the public implementation.

Together they report about **1.6×–7.2×** wall-clock speedup on human / plant / primate whole-genome alignments, without compromising accuracy.

| Input regime | Typical speed gain |
| --- | ---: |
| chromosome-scale assemblies | large (often several-fold) |
| many small contigs | smaller (~1.5× range) |
| ordinary long-read mapping | often modest / comparable to other accelerated minimap2 forks |

Example from the paper’s reporting: barley assembly-to-reference on 48 threads drops from ~12 h (minimap2) to <2 h (mm2-plus).

---

## What did **not** change

Crucially:

* same seeding model (minimap2-style),
* same scoring / chaining logic as the scientific model,
* same intended output semantics / formats,
* presented as near drop-in replacement with near-identical output.

So this is essentially:

```text
same algorithm + better CPU utilization
```

Not a new mapper paradigm. Analogous to:

```text
bwa → bwa-mem2
minimap2 → mm2-plus
```

Nice engineering. Not a methodological revolution.

---

## When it actually helps

You benefit most if you do:

### Assembly-to-assembly alignment

* genome scaffolding helpers (RagTag-style workflows)
* whole-genome synteny
* pangenome / chromosome-scale comparisons
* large complete assemblies against each other

These are exactly the “few long contigs, many threads idle” cases.

### Large genomes + many CPUs

If you routinely run:

```text
minimap2 -t 48
```

on ~chromosome-number query sequences, mm2-plus is aimed at you.

---

## When it probably will not matter

Little reason to rewire pipelines for:

* ordinary ONT / PacBio **read mapping**,
* small genomes,
* highly fragmented assemblies with many contigs,
* jobs already dominated by I/O or downstream analysis rather than chaining.

Classic minimap2 already parallelizes well when there are many independent query sequences.

---

## Pragmatic recommendation

For mixed work that includes repeats, genome comparisons, and assembly-level analysis — but is not dominated by chromosome-vs-chromosome alignment:

1. **Do not switch the whole pipeline.**
2. Keep **minimap2** as the default tool.
3. Optionally install **mm2-plus**.
4. Use it when aligning large, contiguous assemblies.
5. Benchmark once on your real jobs:

```text
/usr/bin/time -v minimap2 ...
/usr/bin/time -v mm2-plus ...
```

If you see clear speedup (say >2×) with equivalent outputs, keep it for those jobs.  
If not, forget it.

Meta-point: this is a common HPC paper pattern — **same algorithm, improved parallelism**. Treat it as a performance patch.

---

## What could actually replace minimap2 later?

Most “new mappers” are also speed tweaks. A few neighboring tools change the problem enough to matter in special cases:

| Tool | Real difference | When it matters |
| --- | --- | --- |
| **Winnowmap** | weighted minimizers; still minimap2-family | highly repetitive genomes, centromeres, satellite-rich assemblies |
| **MashMap3** | MinHash sketch approximate mapping | ultra-fast genome–genome homology / ANI-style scans; not base-level VCF mapping |
| **MMseqs2** | heavy k-mer prefiltering / search engine | huge database search and clustering; not a general long-read mapper |
| **VG / GraphAligner** | pangenome graph alignment | population / graph genomics; currently harder and often slower |

Reality check: minimap2 remains the default because it balances speed, robustness, and flexibility. Specialized tools sit around it. Graph aligners may matter more over the next decade for population genomics, but they are not yet a casual replacement.

---

## Key takeaways

* mm2-plus accelerates minimap2 for chromosome-scale whole-genome alignment mainly by finer-grained parallel chaining and faster chain classification (plus sorting / SIMD engineering).
* Alignment science stays the same; this is not a new mapper class.
* Big wins when few long contigs waste threads; small wins for ordinary read mapping / fragmented assemblies.
* Default advice: don’t switch pipelines; optionally use mm2-plus where assembly-to-assembly wall time hurts.
* If you need a different tool, choose it for a different problem (repeats → Winnowmap; ultra-fast genome compare → MashMap; graphs → VG/GraphAligner), not because a paper said “faster.”
