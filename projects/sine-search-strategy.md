---
title: "SINE search: why fast aligners miss copies, and what to run instead of plain ssearch36"
type: notes
---

> Type: notes (measurements are reproducible with `projects/sine-search-bench/`;
> the recommendations are my reading of them)

Question (2026-10-11): SINEderella / `sear` use ssearch36 because it is the
most sensitive way to find SINE copies (100–500 bp, 10–30 % diverged, 10^4–10^6
copies per genome) in a genome. blastn, BLAT, bwa mem and minimap2 miss copies.
Is there a faster strategy with the same recall, a stepwise scheme (cheap tool
first, ssearch36 on what remains), or a case for writing a new aligner?

Sources: the program digests in `.claude/skills/` (ssearch36-fasta36,
ncbi-blast-core, blat-kent, bwa, minimap2, sear2k, sinederella), the survey in
[`methods/sensitive-search-for-diverged-repeats-landscape.md`](../methods/sensitive-search-for-diverged-repeats-landscape.md),
and the benchmark below.

## 1. Bottom line

1. **No seed-and-extend tool reaches ssearch36 recall above ~25 % divergence,
   and most of them also break on copy number.** The failure is structural
   (seed length, occurrence caps, per-target output caps), not a tuning issue.
2. **Two tools match ssearch36 recall on the test sets: fasta36 with ktup 4,
   at about one seventh of the time, and nhmmer, at the same time as ssearch36
   with a single-sequence query or half of it with a profile built from real
   copies.** Both report every copy and both are drop-in for `sear`'s
   `-m 8C` parsing (nhmmer needs a small parser change).
3. **Stepwise search works for the right reason.** ssearch36's cost in `sear`
   has two terms: the SIMD scan, which is constant per query and genome, and a
   scalar re-alignment of every window that scores a hit. For a family with
   10^5 copies the second term dominates (your own measurement in
   `docs/MASKING.md`: 37–46 s with 1,500 planted copies vs ~4 s with them
   masked). A fast first pass that finds and masks the easy 85–95 % of copies
   therefore makes the ssearch36 residual pass cheap. Section 5 measures this on
   a real genome.
4. **Do not write a new aligner.** The sensitivity you want is full
   Smith-Waterman or a profile HMM, and both already exist in tuned SIMD form.
   What is worth building is the orchestration around them (section 6).

## 2. Why each fast tool misses SINE copies (from the source)

| Tool | Seed | Occurrence cap | Output cap | Consequence for a 357 bp SINE with 10^5 copies at 10–30 % divergence |
|---|---|---|---|---|
| bwa mem | SMEM, exact, ≥19 bp (`-k`) | `-c 500` keeps 500 of n occurrences | secondaries only with `-a`; XA dropped when > 5 hits | Copies below ~85 % identity rarely contain a shared 19-mer; the index is ~5.4 GB for 3 Gbp; scoring `-A1 -B4` has zero expected score at 80 % identity (`bwamem.c:57`). Not a repeat finder at any setting. |
| minimap2 | minimizers k15 w10 | `-f 2e-4` → `mid_occ`, clamped to ≥10; `-f 0` still means 10^6; per streak of repeat seeds only ~1 seed rescued (`seed.c:56-96`) | `-N 5` secondaries, `-p 0.8` | A query made entirely of high-copy seeds keeps one seed; chains need `-n 3` anchors and `-m 40` score. Tuned to k11 w3 -f0 -n1 -m20 -s20 -N100000 it reaches 56 % recall and dies above 25 % divergence. |
| BLAT | 11-mer tiles on a stepSize grid | `-repMatch 1024` drops every tile with ≥1024 occurrences (`genoFind.c:744-769`) | **16 alignments per (query, strand, target sequence)** hard-coded (`gfBlatLib.c:22`) | On a chromosome with thousands of copies you get 32 PSL lines, whatever the options (observed). Default `-minIdentity 90`. fuzzyFind "not good for distant homologies" (`fuzzyFind.h:8-11`). |
| blastn | exact words: megablast 28, blastn 11, blastn-short 7; dc-megablast 11-of-21 template | none for blastn (`db_filter` is magicblast-only) | `-max_target_seqs`, `-max_hsps` | No copy-number problem, but a 300 bp alignment at 70 % identity rarely contains an exact 11-mer; word 7 pushes the cliff from ~25 % to ~30 % divergence. Needs `makeblastdb` or loads `-subject` whole into RAM. |
| LAST | adaptive spaced seeds, suffix array | `-m 10` seed multiplicity default | `-K`/`-C` | Default settings return almost nothing for a high-copy family (59 of 2,789 copies); `-m 1000 -D 1e9` reaches blastn-w7 level; MAM8 seeds cost 100 s of indexing for 47 Mbp and were *less* sensitive here than `-m 1000`. |
| ssearch36 | none (full SW, 8-bit SIMD, 16-bit fallback) | none | one score per 150 kb chunk per strand in the scan; copies appear through recursive sub-alignment (`wm_align.c:149-297`) with threshold `e_cut_r` | Recall limited only by score statistics. Cost = scan + scalar Gotoh over the whole chunk for every hit (`sw_walign`, `wm_align.c:37`). |
| fasta36 | ktup 1–6 lookup, diagonal "init1/initn" regions, banded SW on candidates | none | same as ssearch36 | At ktup 4 a 357 bp query at 70 % identity still has ~30 shared 4-mers on the right diagonal; final alignments are SW. |
| nhmmer | profile HMM, SSV/MSV/Viterbi/Forward filter cascade over windows | none | none | Sensitivity from position-specific scores; the filter cascade is what makes it fast. |

Two `sear`-specific facts from the digests that matter for the strategy:

- `sear:396` passes `-g -3`, which FASTA36 reads as gap-extension −3 (not
  "gapped, both strands"), and a single `-E 2`, which sets the secondary-HSP
  threshold `e_cut_r` to 1.0 (`initfa.c:1126,1138`). With 2.5 kb windows this
  barely matters, but on long chunks it hides weak second copies. Use
  `-E "e1 e2"` explicitly.
- ssearch36 reports at most one score per chunk per strand in the scan;
  `sear`'s 2.5 kb windows with 500 bp overlap are what make every copy
  visible, at the price of 1.25× the genome and boundary duplicates.

## 3. Benchmark design

`projects/sine-search-bench/simulate.py` plants copies of a consensus into a
background sequence and records their coordinates. Copies get a divergence
drawn uniformly from 0–40 % (substitutions with 2:1 transition bias, indels at
12 % of the substitution count, geometric lengths), 20 % are 5′-truncated,
tails of 0–30 extra A. Two consensuses: a random 300-mer (no accidental
homology anywhere) and the real dar_squam1 (357 bp, from
`darevskia-sineome-data`). Background: human chr21 (46.7 Mbp, GRCh38,
unmasked), which has real low-complexity DNA and real repeats but no lizard
SINE. About 2,800 copies per set, at least 1 kb apart.

Real divergence of dar_squam1 copies, measured on the 3,000-copy SubFam sample
(edit distance of the consensus inside the copy): median 10 %, 10th percentile
6 %, 90th percentile 26 %, with a tail to 40 % from old or truncated copies. So
the 20–35 % bins of the benchmark are the ones that decide recall.

Scoring: a planted copy counts as found if a hit covers ≥ 50 % of it. False
positives are hits that overlap no planted copy. Times are wall-clock with 4
threads on this container; index build time is listed separately. ssearch36
and fasta36 search a chunked genome (5 kb windows, 500 bp overlap), the others
the raw genome. Commands are in `bench.py`.

## 4. Results on chr21 + planted copies

Recall per divergence bin (fraction of planted copies found). `fp` = hits
outside planted copies.

### Random 300 bp consensus, 2,799 copies

| tool | index s | search s | hits | fp | recall | 0–5 | 5–10 | 10–15 | 15–20 | 20–25 | 25–30 | 30–35 | 35–40 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ssearch36 | 0 | 8.1 | 3228 | 0 | 0.998 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.99 |
| fasta36 ktup 2 | 0 | 6.2 | 3211 | 0 | 0.994 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.96 |
| **fasta36 ktup 4** | 0 | **1.2** | 3171 | 0 | **0.982** | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.99 | 0.87 |
| nhmmer (single seq) | 0 | 3.3 | 2788 | 0 | 0.996 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.97 |
| nhmmer --max | 0 | 14.4 | 2788 | 0 | 0.996 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.97 |
| blastn word 7 | 0 | 0.9 | 2443 | 0 | 0.822 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.98 | 0.57 | 0.04 |
| blastn word 11 | 0 | 0.6 | 2217 | 0 | 0.765 | 1.00 | 1.00 | 1.00 | 1.00 | 0.98 | 0.81 | 0.32 | 0.02 |
| dc-megablast | 0 | 0.6 | 1660 | 0 | 0.589 | 1.00 | 1.00 | 1.00 | 0.93 | 0.55 | 0.25 | 0.04 | 0.00 |
| megablast | 0.9 | 0.2 | 887 | 0 | 0.317 | 1.00 | 0.91 | 0.45 | 0.12 | 0.03 | 0.00 | 0.00 | 0.00 |
| LAST -m 1000 -D 1e9 | 16 | 0.2 | 2191 | 0 | 0.775 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.90 | 0.31 | 0.00 |
| LAST MAM8 -m 100 | 101 | 0.1 | 1325 | 0 | 0.473 | 1.00 | 1.00 | 0.93 | 0.62 | 0.21 | 0.06 | 0.01 | 0.00 |
| LAST default | 16 | 0.0 | 12 | 0 | 0.004 | 0.03 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| minimap2 k11 w3 -f0 -n1 -m20 | 0 | 6.2 | 1592 | 0 | 0.564 | 1.00 | 1.00 | 0.99 | 0.93 | 0.52 | 0.13 | 0.01 | 0.00 |
| minimap2 -x sr | 0 | 3.0 | 303 | 0 | 0.108 | 0.77 | 0.03 | 0 | 0 | 0 | 0 | 0 | 0 |
| bwa mem -a -k11 -c100000 | 54 | 0.2 | 340 | 0 | 0.121 | 0.85 | 0.06 | 0 | 0 | 0 | 0 | 0 | 0 |
| BLAT default | 0 | 1.3 | 35 | 0 | 0.011 | 0.08 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| BLAT tile 8 step 4 repMatch 10^6 | 0 | 1.6 | 36 | 0 | 0.011 | 0.08 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

### Real dar_squam1 consensus, 2,789 copies

| tool | index s | search s | hits | fp | recall | 0–5 | 5–10 | 10–15 | 15–20 | 20–25 | 25–30 | 30–35 | 35–40 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ssearch36 | 0 | 8.9 | 3266 | 0 | 1.000 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| fasta36 ktup 2 | 0 | 6.6 | 3267 | 0 | 0.997 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.98 |
| **fasta36 ktup 4** | 0 | **1.2** | 3247 | 0 | **0.993** | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.95 |
| nhmmer (single seq) | 0 | 9.3 | 2789 | 0 | 1.000 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| **nhmmer, HMM from 400 real copies** | 0 | **4.3** | 2781 | 0 | 0.996 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.97 |
| blastn word 7 | 0 | 1.2 | 2514 | 0 | 0.837 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.98 | 0.64 | 0.08 |
| blastn word 11 | 0 | 0.8 | 2272 | 0 | 0.784 | 1.00 | 1.00 | 1.00 | 1.00 | 0.98 | 0.89 | 0.41 | 0.03 |
| dc-megablast | 0 | 0.7 | 1680 | 0 | 0.597 | 1.00 | 1.00 | 1.00 | 0.93 | 0.66 | 0.26 | 0.06 | 0.00 |
| LAST -m 1000 -D 1e9 | 20 | 0.2 | 2364 | 0 | 0.830 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.97 | 0.62 | 0.07 |
| LAST MAM8 -m 100 | 96 | 0.2 | 1652 | 0 | 0.590 | 1.00 | 1.00 | 1.00 | 0.94 | 0.63 | 0.25 | 0.05 | 0.00 |
| minimap2 k11 w3 -f0 -n1 -m20 | 0 | 8.7 | 1571 | 0 | 0.562 | 1.00 | 1.00 | 1.00 | 0.94 | 0.58 | 0.12 | 0.02 | 0.00 |
| bwa mem -a -k11 | 52 | 0.5 | 277 | 0 | 0.099 | 0.86 | 0.06 | 0 | 0 | 0 | 0 | 0 | 0 |
| BLAT default | 0 | 1.3 | 37 | 0 | 0.011 | 0.10 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |

Reading the tables:

- Every seed-based tool has a recall cliff whose position is set by its seed:
  megablast ~10 %, dc-megablast and minimap2 ~20 %, blastn word 11 ~27 %,
  blastn word 7 and LAST ~30 %. Real dar_squam1 copies have 10 % of their mass
  beyond 26 % divergence, so even blastn word 7 loses several percent of real
  copies and nearly all of the old ones.
- BLAT and bwa fail for a different reason: copy number. BLAT's 16-per-target
  cap and bwa's seed occurrence handling make them useless for a repeat family
  regardless of divergence.
- fasta36 ktup 4 is the surprise: 7× faster than ssearch36 with 99 % recall
  including the 30–35 % bin, because it only needs a few shared 4-mers on one
  diagonal to trigger a real SW alignment. Its hits are SW alignments with the
  same `-m 8C` columns.
- nhmmer finds everything. With a single-sequence query it costs as much as
  ssearch36; with a profile from 400 real copies it is twice as fast, because
  the profile is more specific and fewer windows survive the filters. The
  profile also rescues old copies that a consensus cannot (your HANDOFF.md
  already records 9 of 145 "nomain" copies rescued this way).
- Zero false positives everywhere at E ≤ 10^-3: with a 47 Mbp background the
  statistics are not the problem; seeding is.

## 5. Real genome: D. valentini (GCA_034642135.1, 1.67 Gbp, 3,559 scaffolds)

*Results pending; filled in below when the run finishes.*

## 6. Recommended strategy

### 6.1 Replace the scan engine inside `sear`, keep everything else

`sear:396` is one line. Three candidates, in order of how little changes:

1. **fasta36 ktup 4 in place of ssearch36** for the genome scan. Same binary
   package, same flags, same `-m 8C` columns, so the awk/bedtools pipeline is
   untouched. Expect ~7× less search time at ≥ 99 % of ssearch36 recall on
   copies up to 35 % divergence. Keep ssearch36 for step 2 (voting on extracted
   copies) and step 3 where libraries are small.
2. **Two-pass: fasta36 ktup 4, mask its hits to N, ssearch36 on the masked
   genome.** The second pass pays the scan term only, since almost no window
   hits, and recovers whatever ktup 4 missed. Total time should be close to one
   ssearch36 scan without copies (section 5 measures it).
3. **nhmmer with a per-family HMM** built from the previous run's copies
   (`hmmbuild --dna` on a MAFFT alignment of a few hundred copies, 1 s). Best
   for old subfamilies and for the "nomain" copies; needs a `--tblout` parser
   in place of the `-m 8C` awk. Run it as an additional pass, not a
   replacement, until its copy sets have been compared on a few genomes.

### 6.2 Hygiene fixes in the ssearch36 call, independent of the above

- Pass `-E "2 1"` (or whatever pair you mean) instead of `-E 2`; the second
  value is the secondary-HSP threshold and the DNA default is 0.001.
- `-g -3` sets gap-extend to −3; if the intent was "both strands", drop it, both
  strands are the default. If the intent was the gap penalty, say so in
  `MANUAL.md`.
- `-Z` is honoured for DNA since July 2025 (`doc/readme.v36:34-38`); pass the
  window count so E-values mean the same for every part and every genome.
- Run consensuses in parallel processes with fewer threads each rather than one
  consensus at a time with `-T nproc`: ssearch36's work unit is a 6 MB buffer of
  windows and a single 300 bp query does not keep 32 threads busy
  (`comp_lib9.c:821-841`).

### 6.3 What not to do

- Do not build a k-mer or minimizer prefilter yourself: pairs at 70–85 %
  identity share no 12–20-mers (your SATELLITES.md measurement), and the
  benchmark shows every such tool loses the 25–35 % bins. fasta36 ktup 4
  already is the right prefilter, with SW alignment attached.
- Do not expect GPU full-SW tools to help today: CUDASW++4.0 is protein and
  score-only, ADEPT returns one alignment per pair, SWIPE needs the genome
  chunked into database records like `sear` already does (survey note,
  finding 6).
- BLAT, bwa mem and minimap2 are not candidates for the scan at any setting.
