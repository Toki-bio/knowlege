---
title: TE-insertion anchor graph — pan-SINEome and Alu topology
type: notes
---

> Type: notes

# One data model for insertion-anchored genomes: pan-SINEome nodes + Alu-graph edges

Working notes (October 2026) connecting three things in this knowledge base and in
[SINE_orth_loc](https://github.com/Toki-bio/SINE_orth_loc):

* the **pan-SINEome** — a reference-free list of SINE insertion loci with the orthologous copy
  or empty site in every genome (Darevskia, 7 assemblies);
* the [Alu-anchor graph spec](https://toki-bio.github.io/knowlege/projects/alu-anchor-graph-topological-profiling/)
  and its [evaluation](https://toki-bio.github.io/knowlege/projects/alu-graph-map-sv-proposal/);
* the [sq2 distant orthologous loci pilot](https://toki-bio.github.io/knowlege/projects/sq2-distant-orthologous-loci/).

Status: the shared model and the edge layer are implemented and run on Darevskia; everything
about cancer / Alu below is reasoning, not a result.

---

## The shared model

The Alu spec writes the genome as \(G = (V, E, P)\). The pan-SINEome is already most of that:

| Alu spec | pan-SINEome (`sine_registry.py build`) |
| --- | --- |
| \(V\): Alu loci with coordinates, subfamily, strand | orthologous groups (stable `PSG` IDs) with the copy or empty site in every genome, family, subfamily |
| \(P(v_i\ \mathrm{present})\) across haplotypes | the group × genome state matrix (P / A / U / X); column frequencies are presence probabilities |
| node presence from short reads | `sine_genotype.py` (genotypes at groups from reads) |
| \(E\): spacers between neighbouring Alus | **new:** `edges.tsv` — groups whose sites are neighbours along each genome, spacer per genome |
| outlier edges | **new:** `breakpoints.tsv` — adjacencies of one genome broken in another; `dupblocks.tsv` |

So the general object is an **anchor graph of TE insertions**: nodes = orthologous insertion
loci (present or empty), edges = adjacency along each genome, states from assemblies or from
reads. The same model serves three timescales:

| Timescale | Question | Main signal |
| --- | --- | --- |
| between species (Darevskia, sq2) | phylogeny, hybrid origin, ancient insertions | node states |
| between individuals / assemblies | polymorphism, misjoins, haplotigs | node states + edges |
| somatic (Alu cancer spec) | rearrangement load | edges |

Share the code and data model; keep them **separate projects** — questions, data and validation differ.

---

## What the SINE work changes in the Alu spec

### 1. Uniqueness belongs to the junction, not the Alu body

The spec defines \(W_u\) on the Alu interval, and the evaluation therefore drops young AluY
as unusable anchors. SINE_orth_loc and `sine_genotype.py` never need the element body to be
unique: they anchor on the **unique flank at the insertion junction**. A read crossing the
junction places uniquely even for near-identical copies (the same principle as MELT-type
insertion callers).

Defined on junction windows (flank k-mers on both sides of the insertion), young polymorphic
AluY become the **most** informative nodes, not the worst. Keep body divergence as a
subfamily/age annotation, not as the anchor filter.

### 2. Short reads are enough for nodes — measured

Real-data test (D. valentini individual 245, the animal behind the dvl assembly; pan-SINEome
built *without* dvl; ~180M Illumina reads, 10x library, read from SRA):

* 92,277 homozygous calls agree with its own assembly, **343 contradict (0.37%)**;
* 9.7k loci unresolved in the assembly get a homozygous call from reads;
* ~1% heterozygous calls (expected: polymorphic insertions, haploid assembly).

The 0.37% is an upper bound on genotyping error — part of it is assembly error at those loci.

### 3. Edges need spanning molecules — linked reads are the cheap middle tier

Short reads cannot measure spacers longer than the insert. Between short reads and ONT there
is a cheaper tier: **linked reads** (10x Chromium — discontinued — or TELL-seq / stLFR): barcodes
span tens of kb, enough to support or reject an adjacency between two anchors without sequencing
across the spacer. The 245 data used above is exactly such a library.

### 4. Multimappers need an explicit policy

The SINE pipeline's "ambiguous:2" class (a copy whose flanks place equally well in two loci)
and the two-flank rescue are the practical form of the spec's low-\(W_u\) and DI/ID nodes. The
spec should say what happens to a node with two equal placements: drop, keep as a fork, or
resolve with spanning reads.

---

## Critique beyond the existing evaluation

1. **Spacer Z cannot see copy number.** A whole-arm gain leaves every spacer the same length.
   \(H_\text{norm}\) counts breakpoints inside edges (≤ 50 kb), not CIN in general.
2. **The real competitor is low-pass copy-number CIN.** Shallow WGS copy-number scores
   (ichorCNA-type, genomic-scar scores) are cheap, established and work on cfDNA. The COLO829
   PoC should report tumor-vs-normal separation **against that baseline**, not only against the
   normal.
3. **COLO829 does not need the pangenome baseline.** With a matched normal, tumor edges can be
   compared with the normal's. Test "matched tumor/normal" and "unmatched vs pangenome" as two
   separate claims.

---

## First implementation and what it shows (Darevskia)

`sine_registry.py build` now writes `edges.tsv`, `breakpoints.tsv` (with `a_specific` = adjacency
kept by no other genome) and `dupblocks.tsv` (runs of ≥ 3 multicopy sites). Simulation: with no
rearrangements every adjacency is kept; a faked misjoin gives exactly one genome-specific
breakpoint plus the two broken flanking adjacencies.

On the 7 Darevskia assemblies (dar7: 149,234 groups, 247,529 edges):

* mix, unp and unm keep 97–98% of each other's adjacencies; dva–mix 96%;
* **arm and nai** each carry 7–9k genome-specific local order changes (sites 10 kb–1 Mb apart
  elsewhere become neighbours) — candidate scaffolding errors (arm is a collapsed assembly of a
  hybrid: haplotype switching is one candidate cause) or real rearrangements;
* **dva** has 40 duplicated blocks (131 sites) vs 1–4 in most others — consistent with regions
  assembled twice (independently: 6,997 mix copies place equally well in two dva loci).

Data: [darevskia-sineome-data `pan/`](https://github.com/Toki-bio/darevskia-sineome-data/tree/main/pan).

This is the cheap test of the graph machinery before any human / cancer use: if it finds
real misjoins and haplotigs here (check against Hi-C maps / read support), porting to HPRC or
COLO829 is a change of inputs, not of method.

---

## Assembly QC from reads at SINE loci

Reads from the **same individual** as an assembly turn the anchor set into ~100k checkpoints:

* reads contradict the assembly (P called 0/0 or A called 1/1) → misassembly, collapsed repeat,
  dropped copy, or genotyper error — inspect gaps, contig ends, depth, spanning pairs at the locus;
* duplicated blocks: true duplication → full depth at both copies; haplotig → about half each;
* heterozygous calls mark where phased assemblies must differ.

It localises and confirms errors; fixing them is polishing / reassembly (Merqury, Inspector,
Flagger, purge_dups, Hi-C remain the genome-wide tools).

---

## Bottom line

```text
nodes   = orthologous insertion loci (pan-SINEome groups), present or empty
edges   = adjacency along each genome, with spacer
anchor  = the junction flank, not the element body (young copies are good nodes)
reads   = nodes from short reads (0.37% contradiction on real data); edges need spanning molecules
first   = Darevskia: misjoins, haplotigs, rearrangements — then HPRC / COLO829
```
