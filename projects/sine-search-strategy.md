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

1. **No seed-and-extend tool reaches ssearch36 recall on old copies, and most
   also break on copy number.** The failure is structural (seed length,
   occurrence caps, per-target output caps), not a tuning issue. On real
   dar_squam1 copies, which are mostly 85–95 % identical, blastn word 7 still
   recovers 99.3 % of ssearch36's loci, but loses 22 % of the 60–70 % ones.
2. **fasta36 with ktup 4 is a drop-in replacement for the genome scan: 99.8 %
   of ssearch36's copies on the real genome, identical output columns, 8×
   faster (51 s vs 403 s for 1.67 Gbp on 4 threads).** The 0.2 % it misses
   are copies below 70 % identity.
3. **A profile HMM (nhmmer) is the only thing that finds copies ssearch36
   misses**: 416 extra loci on D. valentini and 94 % of the 60–70 % bin, at
   half the ssearch36 time. Worth running as an additional pass for old
   subfamilies, not as the primary scan.
4. **ssearch36 has a hard cap of 60,000 retained scan hits per run**
   (`MAX_BEST`, `defs.h:104`); when it fills, the weakest quarter is dropped.
   This is the mechanism behind the hit loss you saw in flankscan. One
   ssearch36 call over a whole genome keeps 43 % of the copies. `sear`'s 20
   parts are what protect SINEderella from it; keep them.
5. **Stepwise search works, for a reason you measured yourself.** ssearch36's
   time has a scan term and a per-hit scalar re-alignment term; for a 10^5-copy
   family the second dominates. A fast first pass that finds and masks the easy
   copies makes the ssearch36 residual pass cheap (section 5).
6. **Do not write a new aligner.** The sensitivity you want is full
   Smith-Waterman or a profile HMM, and both exist in tuned SIMD form. What is
   worth building is the orchestration around them (section 6).

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

Query dar_squam1, 4 threads, same commands as above. No ground truth, so the
reference is part-wise ssearch36 (20,000-window parts, as `sear` does) and
"recall" is the fraction of its merged hits that another tool also covers
(≥ 50 % overlap; HSPs within 50 bp merged first). Identity bins are
ssearch36's gapped percent identity of the reference hit.

| tool | search s | merged hits | recall vs ssearch36 | 95–100 | 90–95 | 85–90 | 80–85 | 75–80 | 70–75 | 60–70 |
|---|---|---|---|---|---|---|---|---|---|---|
| ssearch36, 19 parts | 403 | 141,026 | ref | 17,549 | 75,123 | 30,150 | 12,115 | 15,867 | 9,759 | 2,387 |
| **fasta36 ktup 4, 19 parts** | **51** | 141,224 | **0.998** | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 0.99 | 0.93 |
| nhmmer, HMM from 400 real copies | 179 | 140,924 | 0.994 (+416 not in ref) | 0.98 | 1.00 | 0.99 | 1.00 | 1.00 | 0.99 | 0.94 |
| nhmmer, consensus as query | 170 | 138,595 | 0.981 | 0.97 | 0.99 | 0.98 | 0.97 | 0.97 | 0.96 | 0.88 |
| blastn word 7 | 31 | 141,275 | 0.993 (+283) | 1.00 | 1.00 | 1.00 | 1.00 | 0.99 | 0.97 | 0.78 |
| dc-megablast | 28 | 134,871 | 0.950 | 0.99 | 0.99 | 0.98 | 0.92 | 0.91 | 0.73 | 0.36 |
| megablast | 16 (+16 index) | 116,360 | 0.812 | 0.99 | 0.98 | 0.90 | 0.69 | 0.35 | 0.11 | 0.05 |
| LAST -m 1000 -D 1e9 | 2 (+819 index) | 42,346 | 0.296 | 0.65 | 0.40 | 0.21 | 0.08 | 0.02 | 0.01 | 0.01 |
| ssearch36, **one library of 370k windows** | 190 | 61,435 | 0.431 | 0.79 | 0.68 | 0.08 | 0.08 | 0.09 | 0.08 | 0.07 |
| fasta36 ktup 4, one library | 39 | 60,936 | 0.427 | 0.80 | 0.68 | 0.07 | 0.08 | 0.09 | 0.08 | 0.07 |

What this says:

- **The 60,000-hit cap is real and it is the mechanism behind "ssearch36
  loses hits when the library is huge and one family dominates"** (your
  flankscan HANDOFF.md). `bestp_arr` holds `MAX_BEST` = 60,000 scan hits
  (`defs.h:104`); when it fills, the weakest quarter is dropped and the cutoff
  raised (`compacc2.c:2592-2602`). A single run over 370k windows kept 52,002
  windows with a bit-score floor of 49.5 and lost 57 % of the copies, almost
  all of them below 90 % identity. `sear`'s 20 parts and flankscan's 20k-record
  chunks avoid it by construction; any future "run the whole genome in one
  ssearch36 call" would silently reintroduce it. Keep parts ≤ ~40k windows for
  families of this size, or fewer for bigger ones.
- **Real copies are easier than the planted set.** 66 % of ssearch36's copies
  are above 90 % identity and only 9 % below 75 %, so blastn word 7 already
  recovers 99.3 % of what ssearch36 finds, in 31 s versus 403 s. Its losses are
  concentrated where they hurt the most for subfamily work: 22 % of the 60–70 %
  copies and 3 % of the 70–75 % ones.
- **fasta36 ktup 4 is a drop-in: 99.8 % of ssearch36's copies, same output,
  8× faster.** The 0.2 % it misses are in the 60–70 % bin (old copies).
- **The profile HMM finds copies ssearch36 does not** (416 extra loci) and
  covers 94 % of the 60–70 % bin, at half the ssearch36 time. With a
  consensus as query nhmmer is slightly less complete than ssearch36.
- Index-heavy tools lose on a single query: LAST spends 14 min indexing to
  answer in 2 s; megablast's index is cheap but its seeds are not.

### Two-pass timing (fasta36 ktup 4, mask, ssearch36 on the residual)

*Pending.*

## 6. Recommended strategy

### 6.1 Replace the scan engine inside `sear`, keep everything else

`sear:396` is one line. Three candidates, in order of how little changes:

1. **fasta36 ktup 4 in place of ssearch36** for the genome scan: change the
   binary name and the trailing positional from `3` to `4`. Same package, same
   flags, same `-m 8C` columns, so the awk/bedtools pipeline is untouched.
   Measured: 8× less search time at 99.8 % of ssearch36 recall on D. valentini,
   99.3 % on the planted chr21 set up to 35 % divergence. Keep ssearch36 for
   step 2 (voting on extracted copies) and step 3, where libraries are small
   and the cost is already low.
2. **Two-pass: fasta36 ktup 4, mask its hits to N, ssearch36 on the masked
   genome.** The second pass pays the scan term only and recovers whatever
   ktup 4 missed. Section 5 gives the measured total.
3. **nhmmer with a per-family HMM** built from the previous run's copies
   (`hmmbuild --dna` on a MAFFT alignment of a few hundred oriented copies, one
   second). Best for old subfamilies and the "nomain" copies; needs a
   `--tblout` parser in place of the `-m 8C` awk. Run it as an additional pass
   whose new loci are flagged as HMM-only, not as a replacement, until its copy
   sets have been compared on a few genomes.

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
