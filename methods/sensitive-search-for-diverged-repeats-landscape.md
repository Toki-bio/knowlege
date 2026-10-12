---
title: Tools for sensitive search of diverged repeats, 2013–2026
type: explanation
---

> Type: explanation

Survey of what exists for finding short (100–500 bp), diverged (10–30 %)
repeat copies such as SINEs in gigabase genomes, compiled 2026-10-11 from the
primary sources linked below. Numbers marked *extrapolated* are mine, not the
authors'. See the companion note
[SINE search strategy benchmark](../projects/sine-search-strategy.md) for
measurements on our own test set.

## The landscape in one table

| Tool | Algorithm class | Divergence handled well | Speed class (3 Gbp, one 300 bp query) | All local hits? | Index memory |
|---|---|---|---|---|---|
| ssearch36 | Full Smith-Waterman, striped SIMD, one query vs each library record | Any; limit is score statistics, not seeding | Slowest; no published DNA GCUPS | One best alignment per library record by default (so chunk the genome) | None |
| cross_match -s | Banded SW after word seeding | ~30–35 % | ~30× slower than rmblastn | Yes, masklevel-filtered | Low |
| rmblastn | BLAST seed-extend with cross_match scoring and masklevel | ~25–30 % | Fastest RepeatMasker engine; ~12 h per 3 Gbp on 16 cores for the whole Dfam library (*extrapolated*) | Yes (HSPs) | ~GB |
| nhmmer / pyhmmer | Profile HMM, SSV → Viterbi → Forward filter pipeline over windows | Best measured sensitivity on old copies when the profile is built from an MSA | ~10× slower than rmblastn per library; one 338-position model vs whole human genome took 12 min single-threaded in 2013 | Yes, all envelopes, both strands | Small; optional FM-index (makehmmerdb) ≈ genome size |
| LAST (lastdb -uMAM8 / MAM4) | Adaptive, transition-constrained spaced seeds on a suffix array, gapped extension | Designed for ~70 % identity (human–mouse) | Minutes to hours per mammal genome | Yes; -D sets E per query letter, -K/-C control overlapping hits | MAM8 ~50 GB, MAM4 ~25 GB for a mammal genome |
| LASTZ | One 12-of-19 spaced seed + transition, HSP threshold | ~70 % | Slower and less sensitive than LAST multi-seed (Frith & Noé 2014) | Yes | Moderate |
| MMseqs2 --search-type 3 | Exact spaced 15-mer prefilter, ungapped, then gapped SW | ≥80 % only | Very fast | Single-diagonal bands; adjacent HSPs not merged | ~8 GB k-mer table + DB |
| minimap2 -x asm20 | Minimizers k19 w10, chaining | ≤10–20 % | Very fast | Chains; secondaries capped | ~GB |
| SWIPE -p 0 | Full SW, inter-sequence SIMD over a BLAST db | Any | 106 GCUPS on 12 cores (protein, 2011); DNA not published | Alignments; genome must be chunked into records | None |
| CUDASW++4.0 (2024) | Full SW on GPU | Any | 1.9–5.7 TCUPS | No: protein only, score-only, one score per record | GPU RAM |
| ADEPT (2020) | Full SW on GPU, DNA and protein | Any | 497 GCUPS DNA on 8 GPUs | One alignment per sequence pair | GPU RAM |
| WFA2-lib, Edlib, block aligner, ksw2 | Pairwise global / ends-free | Refinement only | n/a | No, pairwise | n/a |
| nail (2024) | Sparse Forward/Backward after MMseqs2 seeds | HMM-like | 26× faster than HMMER3 | Protein only today | — |
| HAVAC FPGA (2024) | nhmmer SSV filter on an Alveo U50 | Same as nhmmer | ~65× end-to-end vs 1-thread nhmmer | Filter stage only | FPGA |

## Findings

1. **Profile HMMs beat consensus search on diverged TEs, repeatedly.**
   Dfam 2013, human 516 Mb benchmark: nhmmer covered 53.9 % at 0.06 % FDR,
   cross_match 51.0 % at 0.11 %, rmblastn 49.8 % at 0.08 %, sensitive blastn
   44.8 %. RepeatMasker's current help file (Dfam 3.9, 1 Mbp, 16 cores) reports
   human 49.3 % masked with HMMER vs 43.6 % RMBlast vs 44.6 % cross_match -s;
   mouse 62.3 vs 54.6 vs 54.5; *C. elegans* 34.6 vs 23.2 vs 23.5, with
   reversed-sequence false-positive rate ≤0.6 %. Caveat: HMM runs use Dfam
   profile HMMs while BLAST runs use consensus libraries, so library content
   is confounded with algorithm.
   [Dfam 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3531169/),
   [repeatmasker.help](https://raw.githubusercontent.com/Dfam-consortium/RepeatMasker/master/repeatmasker.help),
   [RepeatMasker issue 330](https://github.com/Dfam-consortium/RepeatMasker/issues/330).
2. **nhmmer is slow but tractable and embarrassingly parallel.** Runtime
   scales with number × length × complexity of profiles, so a handful of SINE
   families is cheap even where the full Dfam library is not. makehmmerdb's
   FM-index mode is ~10× faster "with small loss of sensitivity" and flagged
   experimental. pyhmmer exposes nhmmer with the long-targets pipeline and a
   `cpus` argument.
   [nhmmer paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3777106),
   [makehmmerdb](https://manpages.org/makehmmerdb),
   [pyhmmer nhmmer](https://pyhmmer.readthedocs.io/en/latest/api/hmmer/seq.html).
3. **Contiguous exact seeds are the weakest option at 70–80 % identity.**
   PatternHunter: a weight-11 contiguous seed hits a 64 bp region at 70 %
   identity 30 % of the time; the spaced seed `111*1**1*1**11*111` hits 46.6 %
   at the same random-hit rate. Frith & Noé 2014: exact-match seeds are
   "significantly worse"; transition-constrained multi-seed sets (8 seeds)
   beat LASTZ defaults at under half the runtime; memory 8/14/25/48 GB for
   1/2/4/8 seeds on the dog genome. These became LAST's MAM4/MAM8 schemes.
   [Frith & Noé 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC3985675/),
   [LAST seeds doc](https://gitlab.com/mcfrith/last).
4. **LAST is the strongest practical seed-and-extend option for 20–30 %
   divergence.** `lastdb -uMAM8` (or MAM4 at half the memory) with `-c`;
   `lastal -m` raises seed multiplicity (more sensitive), `-D` sets E-value
   per query letter, `-K`/`-C` keep overlapping hits, `-P` threads. No
   published SINE-specific comparison against ssearch exists.
5. **MMseqs2 nucleotide mode is not a diverged-DNA tool.** Default exact
   spaced 15-mers, adjacent HSPs not merged, `-s` changes k-mer count rather
   than allowing substitutions for DNA.
   [issue 253](https://github.com/soedinglab/MMseqs2/issues/253),
   [issue 1116](https://github.com/soedinglab/MMseqs2/issues/1116).
6. **Accelerated full-SW is protein- or pairwise-oriented.** No 2023–2026
   GPU tool was found that does DNA query vs genome with all local hits
   reported. SWIPE comes closest but needs the genome chunked into database
   records.
   [CUDASW++4.0](https://pmc.ncbi.nlm.nih.gov/articles/PMC11531700/),
   [ADEPT](https://pmc.ncbi.nlm.nih.gov/articles/PMC7493400/),
   [SWIPE](https://github.com/torognes/swipe),
   [WFA2-lib](https://github.com/smarco/WFA2-lib).
7. **TE annotation pipelines all collect copies with rmblastn/BLASTN**:
   RepeatModeler2, Earl Grey (BLAST-Extract-Extend), EDTA, TEtrimmer. None
   uses nhmmer or full SW for copy collection, so their copy sets under-sample
   old copies. TE_Bench (Mobile DNA 2026) compares pipelines but gives no
   per-divergence sensitivity.
   [Earl Grey](https://github.com/TobyBaril/EarlGrey),
   [TEtrimmer](https://pmc.ncbi.nlm.nih.gov/articles/PMC12462492/),
   [TE_Bench](https://link.springer.com/article/10.1186/s13100-026-00405-z).
8. **RepeatMasker speed gears**: `-s` is 0–5 % more sensitive at 2.5× the
   time, `-q` 5–10 % less sensitive at 3–4× faster, `-qq` misses 10–25 %.
   The settings "have less effect on NHMMER".

## What this implies for SINE search

Two-stage, with no new aligner needed as the first bet: build a profile HMM
per family from an MSA of ssearch-recovered copies (`hmmbuild --dna`), run
nhmmer on all cores; cross-check with LAST MAM8/MAM4 at permissive `-m`/`-D`;
refine candidate loci with ssearch36. Estimate the false-discovery rate on a
reversed (not complemented) genome, the Dfam protocol. Whether any of this
actually beats ssearch36 on SINE copies specifically is measured in the
companion benchmark note.
