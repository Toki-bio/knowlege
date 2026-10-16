---
title: "ssearch36-fasta36: Statistics: how E-values are estimated"
type: reference
generated: true
---

> Part of the [ssearch36-fasta36](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/ssearch36-fasta36/references/statistics.md)

{% raw %}
# Statistics: how E-values are estimated

All code in `src/scaleswn.c` (linked as `scale_se.o`); orchestration in
`src/comp_lib9.c:1043-1216`; sample collection in `src/compacc2.c:2496-2545`.

## The sample

- One value per library entry (chunk): the best score over strands
  (`compacc2.c:2490-2545`, saved when `frame == nitt1`). Up to `MAX_STATS` =
  60,000 scores (`defs.h:106`); with `-DSAMP_STATS` (set in
  `make/Makefile.linux64_sse2:32`) later entries are sub-sampled via
  `samp_stats_idx()` (`compacc2.c:1331+`).
- `ldb.entries` counts scored chunks (continuations of a long sequence are
  *not* counted again at read time, `comp_lib9.c:2616-2620`, but save_best
  counts every returned result, `compacc2.c:2490-2494`). `notes`: so a
  200 Mbp chromosome contributes ~1,334 entries to the statistics sample.
- When the sample fills early, `process_hist()` is run during the scan and
  later scores are converted to z-scores immediately (`compacc2.c:2500-2516`).

## `-z` modes (`scaleswn.c:30-40`, dispatch `process_hist()` 236-369)

| `-z` | Name | Method | Fallbacks |
|---|---|---|---|
| -1 | none | no statistics; hits sorted by raw score; `nshow = min(-b, nbest)` (`comp_lib9.c:1384-1385`) | |
| 0 | AVE_STATS | mean/variance only, no length correction (`proc_hist_n` 1576) | |
| **1 (default)** | REG_STATS | weighted regression of mean score vs ln(n1) in length bins, outliers above a Karlin-Altschul derived `max_hscore` are excluded (`calc_thresh` 375-427, `proc_hist_r` 429-537, `fit_llen`) | if `zstrim < 0` (too few scores to trim) -> `-z 3` (449-451); if score range < 10 -> `-z 0` (460-463); if fit fails -> `-z 2` (468-470) |
| 2 | MLE_STATS | maximum-likelihood lambda/K, censoring 250 low/high scores (`proc_hist_ml` 837; `fasta36.1:403-405`) | |
| 3 | AG_STATS | Altschul-Gish precomputed lambda/K for protein matrices (`proc_hist_a` 625, `ag_parm` 712) | protein only in spirit; for DNA uses the generic path |
| 4 | REGI_STATS | regression with iterative outlier removal (`proc_hist_ri` 541) | |
| 5 | REG2_STATS | regression with length-scaled variance (`proc_hist_r` 348) | |
| 6 | MLE2_STATS | MLE with per-sequence composition (`proc_hist_ml2` 903); needs `do_karlin()` per chunk (`dropgsw2.c:832-838`) | |
| 11..16 | shuffled library | every chunk is also shuffled and scored (`buf_shuf_work`, `work_thr2.c:357-359`); the shuffled scores are the sample (`compacc2.c:2531-2534`); method = mode-10 (`scaleswn.c:314-316`) | doubles scan time |
| 21..26 | shuffled top hits | normal E() plus `E2()` from shuffles of the hits with `E < e_cut` (at least 10, each >= 5 times; `comp_lib9.c:1174-1216`); second method from `-z "21 2"` (`initfa.c:1312-1325`) | `zsflag2 == 3` is replaced by 2 (`scaleswn.c:310`) |

Forced overrides inside `process_hist()`:
- `nstats < 10` -> `-z 3` (`scaleswn.c:289`).
- query length `n0 <= 10` with regression modes -> `-z 2` (`291-293`).
- `zdb_size` defaults to `db.entries` (`319`; also `comp_lib9.c:1103`).

## From z-score to E()

- `find_z*()` converts a raw score to a length-normalised z-score (mean 50,
  sd 10) using the fitted `rho`, `mu`, `mean_var` (`find_zr` 2480).
- `zs_to_E(zs, n1, dnaseq, entries, db)` (`scaleswn.c:2681-2720`):
  `z = (zs-50)/10`, extreme-value tail `e = exp(-(z*pi/sqrt6 + Euler))`,
  `E = k * (1 - exp(-e))`. **k** = `entries` (= `zdb_size`) for protein, but for
  DNA/RNA without `-Z`: `k = db.length / n1` (`2694-2699`), i.e. total library
  residues divided by the length of *this* chunk, so with 150 kb chunks k is
  ~the number of chunks. `-Z` for DNA:DNA is honoured only since July 2025
  (`doc/readme.v36:34-38`; `db.zdb_size_set` wired at `comp_lib9.c:567,891`).
- `scale_scores()` assigns `rst.escore` to every kept hit and sorts
  (`2942-2967`); bit scores via `zs_to_bit()` (2668).
- `E1_to_s()` (2787-2830) inverts this to get the raw score for a target E; used
  for `repeat_thresh` (see `library-chunking-and-hsps.md`).

## How `-E` and `-z` interact

- `-E e1` sets `e_cut` (`initfa.c:1113-1141`); after `scale_scores()` the
  displayed set is `escore < e_cut` (`comp_lib9.c:1364-1367`), modified by `-b`
  (`1369-1385`: `-b N` caps, `-b =N` forces exactly N ignoring E, `-b >N`
  floor) and `-F e_low` skips the most significant (`1340-1343`).
- E-values are only as good as the fit; `-z` picks the fit. With `-z -1`
  there are no E-values and `-E` is meaningless (`1384-1385`).
- `-z 2x` adds a second column `E2()`; `-E` still filters on the first.
- `-Xb` shows z-scores instead of bits in the hit list (`initfa.c:1371`).

## Small library: one genome of ~30 chromosomes

Because of chunking the statistics never see "30 sequences": a 3 Gbp genome
gives ~20,000 150-kb chunks, every one scored on both strands, so `nstats` is
~20,000 (<= 60,000) and `-z 1` regression has plenty of points, but **all at
nearly the same length** (149,850-150,000). `notes`: the regression on ln(n1)
then has almost no leverage and degenerates to a mean/variance fit at
n1 = 150 kb; the 150 kb chunk distribution is what the EVD is fitted to. That is
statistically fine for E-values of other 150 kb chunks (which is what is
reported), but the last chunk of each chromosome is shorter and sits off the
fitted line.

Automatic safety nets:
- If fewer than `shuff_max` (500, `-k`, `initfa.c:969,1206-1209`) hits were
  kept and `-z` is not 3 or -1, the program shuffles every kept hit and refits
  with `-z 1x` (`comp_lib9.c:1079-1096`). A genome smaller than ~75 Mbp (<500
  chunks) therefore silently switches to shuffle statistics.
- Fewer than 10 scores -> `-z 3` (`scaleswn.c:289`).

Recommendations (`notes`, derived from the code above and the guide):
- Default `-z 1` is acceptable for a whole vertebrate genome; check the fitted
  line in the header (`Expectation_n fit: rho ... mu ... mean_var`,
  `scaleswn.c:531-534`).
- For a repeat-rich query where many chunks contain true hits (a SINE with
  10^5 copies: *every* chunk has a real match), the "random" distribution is
  contaminated; `-z 11` (score shuffled copies of every chunk) gives a clean
  null at 2x cost, `-z 21` is cheaper and gives a second `E2()` from shuffles of
  the high scorers only (`doc/fasta_guide.tex:1503-1513`; `fasta36.1:411-427`).
- `-Z N` forces the database size used in E() (`doinit.c:488-490`); set it to
  the number of chunks (or genome_length/150000) to make E-values comparable
  across genomes.
- `-E` for a genome scan: the default 2.0 with ~20,000 x 2 strands tested will
  show ~2 false hits; raise sensitivity with `-E 10 -z 11` and filter later.
{% endraw %}
