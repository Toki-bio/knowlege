---
title: ERAST and orthologous SINE loci
type: notes
---

> Type: notes

# ERAST for orthologous SINE loci — verdict: not a fit

[Scalable homology detection with ERAST](https://www.nature.com/articles/s41587-026-03051-1)
(Nature Biotechnology 2026, Zhejiang University / Tencent;
[code](https://github.com/TencentAILabHealthcare/ERAST)). Assessed from the code and README,
for the question "can it speed up multi-genome comparison of orthologous SINE loci"
([SINE_orth_loc](https://github.com/Toki-bio/SINE_orth_loc)).

## What it is

* Embeds each sequence with a language model (DNA: a Caduceus model trained to tell the phylum
  of genomic fragments; protein: ESM2 + Pfam), searches ~10⁹ stored embeddings by similarity,
  re-ranks the best hits. Built for **remote homology**, reported ~50× faster than Foldseek.
* Output: a ranked list of target IDs per query — **no alignments, no identity, no coordinates**.
* PyTorch + model weights, realistically a GPU; PolyForm Noncommercial licence (academic use allowed).

## Why it does not fit orthologous loci

1. **Opposite task.** Orthology needs the one exact locus among near-identical paralogs (the
   multimapper problem). Embedding search pulls similar sequences together — useful for families,
   wrong for picking one locus among near-copies.
2. **Missing outputs.** The pipeline needs the insertion site, flank identity and SINE presence;
   ERAST gives neither coordinates nor alignments.
3. **Wrong bottleneck.** Flank mapping (bwa / minimap2) takes minutes; the time goes to per-cluster
   MAFFT + ComPair.sh, which ERAST does not touch.

## Where it might help (untested)

Assigning very diverged or truncated copies to SINE families, or finding related families
across distant reptile groups. Try nhmmer with consensus profiles or MMseqs2 clustering first —
fast, CPU-only, real alignments and scores.

## What actually makes many genomes cheap

Map each new genome once against a **shared index of representative flanks** of all
pan-SINEome groups instead of all pairs: ~190 pairwise runs for 20 genomes become 20 runs.
See [TE-insertion anchor graph](https://toki-bio.github.io/knowlege/projects/te-anchor-graph-pansineome/).
