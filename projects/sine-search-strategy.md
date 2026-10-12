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
   time is 50 s of SIMD scan plus 350 s of scalar re-alignment of the 141k hit
   windows on D. valentini. fasta36 ktup 4, mask, then ssearch36 on the residual
   takes about 3 min in total, keeps 99.9 % of full ssearch36's loci and adds a
   few, versus 403 s for ssearch36 alone (section 5).
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

| step | time | result |
|---|---|---|
| fasta36 ktup 4, 19 parts | 51 s | 141,224 loci (193,467 raw HSPs) |
| mask HSPs ± 20 bp to N (`mask_hits.py`, pure Python) | ~90 s | 57.0 Mbp masked, 3.46 % of the genome |
| ssearch36, 19 parts, on the masked genome | **50 s** | 155 merged loci: 141 that full ssearch36 also had, 14 new |
| total | ~3 min | 99.9 % of full ssearch36's loci plus 14 |
| full ssearch36, 19 parts, for comparison | 403 s | 141,026 loci |

The residual pass is the cleanest measurement of ssearch36's cost structure:
the SIMD scan of 1.67 Gbp × 357 bp × 2 strands is **50 s on 4 threads** (about
6 GCUPS), and the other **350 s of the full run is the scalar Gotoh
re-alignment of every hit window** (`sw_walign`, `wm_align.c:37`): 141k hits ×
5 kb × 357 bp ≈ 2.5 × 10^11 scalar cells. In `sear` the windows are 2.5 kb,
so that term halves but still dominates for any family with 10^5 copies. This
is what your MASKING.md measured ("the time is the known copies being
aligned") and it is why masking helps so much. fasta36 pays far less per hit
because its DNA alignments are banded (width 16, `initfa.c:1699`).

## 7. Parameter sweep on a 104 Mbp slice of D. valentini (2026-10-12)

Nine largest scaffolds (104 Mbp), reference = the full-genome part-wise
ssearch36 loci falling in those scaffolds (8,184 merged loci). 4 threads;
windows 5 kb / 500 bp overlap unless stated. Recall = fraction of reference
loci covered; last column = recall in the 60–70 % identity bin, where the
differences are.

| variant | search s | merged loci | recall | 60–70 % bin | note |
|---|---|---|---|---|---|
| ssearch36 `-E 1e-3 -z 3` | 24.6 | 8,185 | 1.000 | 1.00 | reference settings |
| ssearch36 `-E 2 -z 11` | 28.1 | 8,263 | 1.000 (+79) | 1.00 | SINEderella's statistics; `-z 11` shuffles cost 14 % |
| ssearch36 `-E 2 -z 3` | 24.0 | 8,370 | 1.000 (+187) | 1.00 | the +187 are weak fragments, see §7.1 |
| fasta36 ktup 1 | 62.2 | 8,203 | 0.998 | 0.91 | degenerate: every position seeds |
| fasta36 ktup 2 | 17.1 | 8,199 | 0.997 | 0.90 | |
| fasta36 ktup 3 | 6.0 | 8,204 | 0.998 | 0.91 | |
| **fasta36 ktup 4** | **3.1** | 8,204 | 0.998 | 0.91 | |
| fasta36 ktup 5 | 1.9 | 8,202 | 0.997 | 0.90 | |
| **fasta36 ktup 6 (DNA default)** | **1.4** | 8,202 | 0.997 | 0.90 | 17× faster than ssearch36 |
| fasta36 ktup 4 `-E 2 -z 11` | 4.9 | 8,284 | 0.998 (+79) | 0.91 | sear-style statistics |
| nhmmer, consensus | 10.1 | 8,102 | 0.988 | 0.90 | |
| nhmmer, HMM from 400 copies | 10.4 | 8,225 | 0.997 (+39) | 0.93 | |
| blastn word 7, 2/−3 | 1.4 | 8,246 | 0.994 | 0.79 | |
| blastn word 9 | 1.2 | 8,240 | 0.994 | 0.79 | |
| blastn word 11 | 1.3 | 8,225 | 0.992 | 0.78 | |
| **blastn word 7, reward 1 penalty −1, gap 2/1** | 1.6 | 8,202 | **0.997** | **0.87** | scoring matters more than word size |
| blastn word 11, 1/−1, 2/1 | 1.5 | 8,187 | 0.995 | 0.86 | |
| blastn word 7, E ≤ 10 | 1.3 | 8,402 | 0.996 (+181) | 0.81 | |
| dc-megablast | 1.4 | 7,822 | 0.949 | 0.35 | |
| megablast | 0.7 (+1.0 index) | 6,642 | 0.800 | 0.04 | |
| LAST `-m 1000 -D 1e9` | 0.3 (+32 index) | 7,862 | 0.953 | 0.31 | |
| minimap2 k11 w3 -f0 (genome indexed) | 13.3 | 7,767 | 0.947 | 0.63 | |
| **minimap2, roles swapped** (consensus indexed, windows as reads) k9 w1 | 3.1 | 8,052 | 0.975 | 0.80 | no occurrence filter can bite |
| minimap2, swapped, k11 w2 | 2.7 | 7,986 | 0.969 | 0.73 | |
| bwa mem -k 11 (genome indexed) | 1.7 (+112 index) | 3,421 | 0.416 | 0.00 | |
| **bwa mem, roles swapped**, k 9 or 11 | 10 | 6,975 | 0.789 | 0.01 | seeds fixed, scoring (`-A1 -B4`) still kills < 80 % |
| BLAT default | 6.4 | 423 | 0.035 | 0.00 | 16-per-target cap |
| **BLAT, roles swapped**, tile 8 | 2.2 | 7,767 | 0.971 | 0.76 | cap is now per window |
| BLAT, swapped, tile 11 | 4.4 | 7,477 | 0.936 | 0.43 | |

Window-size variants (same slice):

WINDOW_ROWS

### 7.1 Reading the sweep

- **ktup barely matters for recall, only for time.** Every ktup from 1 to 6
  recovers 99.7–99.8 % of the reference, all losing the same 9–10 % of the
  60–70 % identity copies. The DNA default ktup 6 is the fastest at 1.4 s,
  17× faster than ssearch36 on the same windows. The loss is not in the
  seeding stage but in fasta36's banded (width 16) alignment and its
  `init1`/`initn` region thresholds, which is why smaller words do not buy it
  back.
- **blastn's problem is scoring, not word size.** Word 7, 9 and 11 give the
  same recall; switching from the default 2/−3 to reward 1 penalty −1 with gap
  costs 2/1 lifts the 60–70 % bin from 0.79 to 0.87 and overall recall to
  0.997 in 1.6 s. The 2/−3 scheme is tuned for ~90 % identity; 1/−1 is the
  classic choice for ~75 %.
- **Swapping roles rescues the index-based tools from their occurrence
  filters.** With the 357 bp consensus as the indexed reference and the 2.5–5
  kb windows as reads, no k-mer occurs more than a few times in the index, so
  minimap2's `-f`, bwa's `-c` and BLAT's `repMatch` never fire and BLAT's
  16-alignments-per-target cap applies per window. BLAT goes from 3.5 % to
  97 %, bwa from 42 % to 79 %, minimap2 from 95 % to 97.5 %. They still trail
  fasta36 and blastn because their chaining and scoring were designed for
  reads, not 70 %-identity repeats.
- **The raised cap works.** A rebuild with `MAX_BEST` = 4,000,000 in `defs.h`
  run over the whole genome as one 370k-window library found 140,402 loci,
  99.6 % of the part-wise reference, in 414 s, versus 43 % for the stock
  binary. So the loss is a compile-time limit, not an algorithmic one; parts
  or a rebuilt binary are equivalent fixes.
- **The HMM-only loci** (§5) split into two kinds. 25 of 465 are strong
  (270–430 bits) but are not misses: they are tandem pairs or dimers that
  nhmmer reports as one 700–780 bp envelope while ssearch36 reports two
  copies 60 bp apart, which my 50 bp merge rule did not join. The other 440
  are short (median 64 bp of model) and weak (26–40 bits) fragments at ~75 %
  identity. ssearch36 does align them (532 alignments at E ≤ 10^-3 when the
  loci are searched in isolation) but at genome-part scale their E-values
  land at 0.1–2, above the cutoff. The profile scores the same fragment higher
  because it weights the conserved positions, so the same 64 bp reaches
  significance. That is the whole HMM advantage here: statistical power on
  fragments, not detection of copies ssearch36 cannot align. `-E 2` in
  SINEderella already admits the comparable ssearch36 fragments (+79 loci).

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
   genome.** The second pass pays the scan term only (50 s for 1.67 Gbp) and
   recovers most of what ktup 4 missed. Measured total about 3 min versus
   403 s, with 99.9 % of the loci plus a handful of new ones. Use
   `bedtools maskfasta` instead of the Python masker to cut the middle step to
   seconds. This is the variant to choose if the 0.2 % matters.
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
