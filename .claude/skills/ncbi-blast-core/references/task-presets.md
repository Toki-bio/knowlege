# `-task` presets (exact values at the pinned commit)

Where defined: `CBlastOptionsFactory::CreateTask`
(`src/algo/blast/api/blast_options_handle.cpp:337-382`, task list
`GetTasks` :212-220) -> `CBlastNucleotideOptionsHandle`
(`src/algo/blast/api/blast_nucl_options.cpp`) and `CDiscNucleotideOptionsHandle`
(`src/algo/blast/api/disc_nucl_options.cpp`); numeric macros in
`include/algo/blast/core/blast_options.h:57-160`; C defaults in
`src/algo/blast/core/blast_options.c` (`BLAST_FillLookupTableOptions` :1129,
`BLAST_FillExtensionOptions` :699, `BLAST_FillScoringOptions` :823).
The `blastn` binary defaults to `megablast` (`blastinput/blastn_args.cpp:54-55`).

| Option | megablast (default) | blastn | dc-megablast | blastn-short | rmblastn |
|---|---|---|---|---|---|
| preset | `SetTraditionalMegablastDefaults` (:115-134) | `SetTraditionalBlastnDefaults` (:70-90) | `CDiscNucleotideOptionsHandle` (disc:48-91) | blastn + overrides (handle :353-360) | blastn + reward/penalty 0 (handle :365-369) |
| lookup type | `eMBLookupTable` (:148) | `eNaLookupTable` (:139) | MB (disc:55-57) | as blastn | as blastn |
| word_size | 28 (:149) | 11 (:140) | 11 (disc:62) | 7 (:359) | 11 |
| template | none | none | type 0 coding, length 18 (disc:60-61) | none | none |
| reward/penalty | 1/-2 (:220-221) | 2/-3 (:203-204) | 2/-3 (disc:87-90) | 1/-3 (:356-357) | 0/0 = matrix scoring, `-matrix` required |
| gap open/extend | 0/0 linear (:218-219) | 5/2 (:201-202) | 5/2 | 5/2 | 5/2 |
| prelim ext / traceback | greedy / greedy (:192-193) | DP / DP (:182-183) | DP / DP (disc:82-83) | DP | DP |
| xdrop_ungap (bits) | 20 (`BLAST_UNGAPPED_X_DROPOFF_NUCL`, C default; MB preset sets only window :173) | 20 (:165) | 20 (disc:69) | 20 | 20 |
| xdrop_gap (bits) | 25 greedy (:189) | 30 (:179) | 30 (disc:79) | 30 | 30 |
| xdrop_gap_final (bits) | 100 (:190) | 100 (:180) | 100 (disc:80) | 100 | 100 |
| gap_trigger (bits) | 27 (:191) | 27 (:181) | 27 (disc:81) | 27 | 27 |
| window_size | 0 (:173) | 0 (:166) | 40 (disc:72) | 0 | 0 |
| off_diagonal_range | 0 | 0 (:167) | 0 | 0 | 0 |
| dust / soft masking | on, mask_at_hash (:157-158) | on | on | **off** (`ClearFilterOptions` :360) | on |
| strands | both (:159) | both | both | both | both |
| evalue | 10 (:254) | 10 (:234) | 10 | **1000** (:358) | 10 |
| hitlist_size (API) | 500; CLI `-max_target_seqs` default 100 | 500 | 500 | 500 | 500 |
| min_diag_separation | 6 (:259) | 50 (:239) | 50 | 50 | 50 |
| perc_identity / qcov | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 | 0 / 0 |
| mask_level | 101 (off) | 101 | 101 | 101 | `-mask_level` (`blast_args.cpp:541-547`) |
| extras | `-use_index` allowed (`blast_args.cpp:3388-3391`) | | `SetTraditionalBlastnDefaults` throws (disc:93-97) | | `-complexity_adjust`, `-min_raw_gapped_score`/`-xdrop_ungap` as raw cutoff (`blast_args.cpp:554-571`) |

Also present: `vecscreen` (`SetVecScreenDefaults` :93-112: 1/-5, gap 3/3,
E 700, search space 1.75e12) - excluded from `blastn -task`
(`blastn_args.cpp:58`). `rmblastn` is a separate binary (`rmblastn_args.cpp`)
but the task name is accepted by `blastn` too (`GetTasks` :219).

Order matters: lookup defaults before initial-word defaults "because default
scanning stride depends on the lookup table type" (:82-84).

## CLI flags that override presets (all in `blastinput/blast_args.cpp`)

| Flag | Lines | Target |
|---|---|---|
| `-word_size` (>= 4) | :163-171, :293-306 | `SetWordSize` |
| `-reward`, `-penalty` (<= 0) | :644-662, :676-681 | scoring |
| `-gapopen`, `-gapextend` | :173-183, :258-271 | scoring |
| `-no_greedy` | :664-667, :684-688 | DP prelim + traceback |
| `-template_type`, `-template_length` | :698-768 | disc template (`eMBWordTwoTemplates` for coding_and_optimal) |
| `-xdrop_ungap`, `-xdrop_gap`, `-xdrop_gap_final` | :209-226, :273-283 | bits |
| `-window_size` (0 = one-hit) | :470-498 | two-hit |
| `-off_diagonal_range` | :500-524 | `scan_range` |
| `-ungapped` | :907-933 | `SetGappedMode(false)` |
| `-dust`, `-filtering_db`, `-window_masker_taxid/_db`, `-soft_masking` | :326-468 | query masks |
| `-lcase_masking` | :2016, `cmdline_flags.cpp:115-116` | FASTA lowercase |
| `-evalue`, `-perc_identity`, `-qcov_hsp_perc`, `-max_hsps`, `-searchsp` | :143-150, :185-207, :228-232, :253-255, :300-318 | hit saving / stats |
| `-dbsize` | :2234-2240, :2568-2570 | eff. lengths |
| `-culling_limit`, `-best_hit_*`, `-subject_besthit` | :3288-3346 | HSP writers |
| `-max_target_seqs`, `-num_alignments`, `-num_descriptions` | :2723-2733, :2892-2975 | `hitlist_size` |
| `-min_raw_gapped_score` | :705-708, :741-743 | `cutoff_score` |
| `-use_index`, `-index_name` | :3350-3411 | megablast dbindex |
| `-db_soft_mask`, `-db_hard_mask` | :2518-2522 | subject masking (BLAST DB only) |
| `-subject`, `-subject_loc` | :2524-2555 | in-memory subjects |
