---
title: "ncbi-blast-core: Core data structures"
type: reference
generated: true
---

> Part of the [ncbi-blast-core](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/ncbi-blast-core/references/data-structures.md)

{% raw %}
# Core data structures

All in `include/algo/blast/core/` unless noted. Options = user intent (bits,
E-values); Parameters = derived raw scores.

## `LookupTableOptions` (`blast_options.h:206-224`)

| Field | Meaning | Default for blastn / megablast |
|---|---|---|
| `lut_type` | `ELookupTableType` (:187-201): `eMBLookupTable`, `eSmallNaLookupTable`, `eNaLookupTable`, `eNaHashLookupTable` (16-mer hash, mapping), `eIndexedMBLookupTable`, `eMixedMBLookupTable` | eNaLookupTable / eMBLookupTable (`blast_nucl_options.cpp:139,148`); the real type is re-chosen at build time |
| `word_size` | exact matches required to trigger extension | 11 / 28 (`BLAST_WORDSIZE_NUCL/MEGABLAST` :67-68) |
| `mb_template_length`, `mb_template_type` | discontiguous template (16/18/21; 0 coding,1 optimal,2 both) | 0 (dc-megablast: 18, 0) |
| `stride` | skip words when indexing (unused by presets, 0) | 0 |
| `db_filter`, `max_db_word_count` | drop query words occurring >= N times in DB; **mapping only, word>=16** | FALSE; `MAX_DB_WORD_COUNT_MAPPER` 30 (:174) |

## Lookup tables (`blast_nalookup.h`)

- `BlastSmallNaLookupTable` (:62-80): `word_length`, `lut_word_length`, `scan_step`,
  `backbone_size = 4^lut_word_length`, `Int2 final_backbone` (one hit or
  -(overflow index)), `Int2 overflow`, `masked_locations` (soft-mask check).
  Limit: < 32767 table entries and query offsets < 32768 (`blast_nalookup.c:171-176`).
- `BlastNaLookupTable` (:131-156): same with `Int4` cells (`NA_HITS_PER_CELL 3`
  :109) and a presence vector `pv`.
- `BlastMBLookupTable` (:238-269): `hashtable[4^lut_word_length]` -> first query
  position, `next_pos[]` chain, `pv_array` (bit per group of cells,
  `pv_array_bts`), `discontiguous`, `template_length`, `template_type`,
  `two_templates`/`hashtable2`, `scan_step`, `longest_chain`.
- `EDiscTemplateType` (:221-235) with the 12 templates drawn at :191-215.
- `LookupTableWrap` (`lookup_wrap.h:50-63`): `lut_type`, `lut`, plus function
  pointers `read_indexed_db`, `check_index_oid`, `end_search_indication`
  (dbindex) and `lookup_callback` (used for soft-mask seed checks).

## `BlastInitialWordOptions` / `BlastInitialWordParameters`

- Options (`blast_options.h:296-304`): `gap_trigger` (bits, 27.0 for nucl
  `BLAST_GAP_TRIGGER_NUCL` :140), `window_size` (0 = one-hit; `-window_size`),
  `scan_range` (`-off_diagonal_range`, default 0), `x_dropoff` (bits, 20).
- Parameters (`blast_parameters.h:95-123`): per-context `BlastUngappedCutoffs`
  {`x_dropoff_init`, `x_dropoff`, `cutoff_score`, `reduced_nucl_cutoff_score`},
  `container_type` (`eDiagArray` or `eDiagHash` for concatenated query > 8000,
  `blast_parameters.c:171-173,219-225`), `nucl_score_table[256]`,
  `matrix_only_scoring` (rmblastn), `ungapped_extension` (TRUE except PHI).

## `BlastExtensionOptions` / `Parameters`

- Options: `gap_x_dropoff` (30 DP / 25 greedy bits), `gap_x_dropoff_final`
  (100 bits), `ePrelimGapExt` (`eDynProgScoreOnly` | `eGreedyScoreOnly`),
  `eTbackExt` (`eDynProgTbck` | `eGreedyTbck`) (`blast_options.h:321-345`,
  defaults `blast_options.c:660-696`).
- Parameters (`blast_parameters.h:129-134`): raw `gap_x_dropoff`,
  `gap_x_dropoff_final` = bits*ln2/min gapped lambda (`blast_parameters.c:455-461`).

## `BlastScoringOptions` (`blast_options.h:461-479`)

`reward`, `penalty`, `gap_open`, `gap_extend`, `gapped_calculation`,
`complexity_adjusted_scoring` (rmblastn), `matrix` (protein / rmblastn only).
Reward=penalty=0 is the rmblastn "matrix scoring" signal (`blast_options.c:878-880`).

## `BlastHitSavingOptions` (`blast_options.h:391-454`)

| Field | CLI | Where applied |
|---|---|---|
| `expect_value` | `-evalue` (10) | cutoff raw score + final reap (`blast_hits.c:1976`) |
| `cutoff_score` | `-min_raw_gapped_score` | overrides E-value cutoff (`blast_parameters.c:876-896`) |
| `percent_identity` | `-perc_identity` | traceback only (`blast_hits.c:993-1001`) |
| `hitlist_size` | `-max_target_seqs` / `-num_alignments` | number of *subject sequences* kept per query |
| `hsp_num_max` | API `SetMaxNumHspPerSequence` (0 = INT4_MAX) | `Blast_HSPListSaveHSP` heap |
| `max_hsps_per_subject` | `-max_hsps` | post-traceback trim (`blast_hits.c:2049`) |
| `culling_limit` | `-culling_limit` | culling writer |
| `min_diag_separation` | API (50 blastn, 6 megablast) | interval-tree containment |
| `do_sum_stats` | (`-sum_stats` hidden for blastn) | FALSE for gapped blastn |
| `query_cov_hsp_perc` | `-qcov_hsp_perc` | post-traceback |
| `mask_level` | `-mask_level` (rmblastn, 101 = off) | `Blast_HSPResultsApplyMasklevel` |
| `hsp_filt_opt` | best-hit / culling / `-subject_besthit` | `BlastHSPFilteringOptions` (:371-389) |
| `min_hit_length` | none found | `s_HSPTest` |

`BlastHitSavingParameters` (`blast_parameters.h:166-185`): per-context
`BlastGappedCutoffs{cutoff_score, cutoff_score_max}`, `cutoff_score_min`,
`link_hsp_params` (NULL for gapped blastn), `prelim_evalue`.

## `BLAST_SequenceBlk` (`blast_def.h:242-290`)

`sequence` (ncbi2na packed for subjects; blastna 1 byte/base for query),
`sequence_start` (sentinel byte before), `length`, `oid`, `sequence_nomask`
(query without filtering for identity counts), `compressed_nuc_seq` (query
packed for fast ungapped extension), `lcase_mask` (query mask), `chunk`,
`seq_ranges`/`num_seq_ranges`/`mask_type` (subject soft/hard masking).

## `BlastQueryInfo` / `BlastContextInfo` (`blast_query_info.h:60-91`)

One context per query strand (blastn: 2 per query). `query_offset`,
`query_length`, `eff_searchsp`, `length_adjustment`, `frame`, `is_valid`.
Queries are concatenated into one search (`GetQueryBatchSize`: blastn 100 kb,
dc-megablast 500 kb, megablast 5 Mb; `blast_input_aux.cpp:104-113`).

## HSP containers (`blast_hits.h`)

`BlastHSP` (:126-148): `score`, `num_ident`, `bit_score`, `evalue`, `query`/`subject`
`BlastSeg` {frame, offset, end, gapped_start}, `context`, `gap_info`.
`BlastHSPList` (:153-166): per subject oid; `hsp_max`, `best_evalue`.
`BlastHitList`: per query, at most `hsplist_max` = prelim hitlist size; heap by
best E-value once full (`blast_hits.c:3243-3297`). `BlastHSPResults`: per query.

## Ungapped hit containers (`blast_extend.h`)

`BLAST_DiagTable` (:77-95): `hit_level_array[diag]` {last_hit, flag},
`hit_len_array` (two-hit), size = next pow2 >= qlen + window
(`blast_extend.c:47-76`). `BLAST_DiagHash` (:98-106) for long queries.
`BlastInitHSP`/`BlastInitHitList` (:150-164) hold `q_off, s_off` and
`BlastUngappedData{q_start, s_start, length, score}`.

## `BlastScoreBlk` (`blast_stat.h`)

`matrix` (16x16 blastna incl. ambiguity codes, `blast_stat.c:1060-1130`),
`kbp_std[ctx]` (ungapped K-A), `kbp_gap_std[ctx]` (gapped, from tables),
`kbp`/`kbp_gap` aliases, `reward`, `penalty`, `round_down` (odd scores
rounded down for some reward/penalty pairs), `matrix_only_scoring`,
`complexity_adjusted_scoring`, `gbp` (Gumbel; NULL for blastn).
{% endraw %}
