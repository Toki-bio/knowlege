# Two-hit window, ungapped extension, gapped extension, HSP saving and limits

## 1. Two-hit window and diagonal bookkeeping (`na_ungapped.c:631-800`)

`s_BlastnDiagTableExtendInitialHit` (hash variant :803 for concatenated
queries > 8000 bases, `blast_parameters.c:171-173,219-225`):

1. `two_hits = (window_size > 0)` (:654). Defaults: blastn 0, megablast 0,
   dc-megablast 40 (`blast_options.h:58-61`; `-window_size`, `blast_args.cpp:470-498`;
   `BLAST_GetSuggestedWindowSize` returns without changing blastn
   `blast_options.c:1211-1217`).
2. A hit already inside the explored span of its diagonal is rejected (:670).
3. Two-hit mode: a new hit becomes `hit_ready` only if the previous hit on the
   same diagonal was saved or ended more than `window_size` before (:672); the
   first hit is recorded with its length (:794) and the second hit within the
   window triggers extension. Off-diagonal second hits are accepted within
   `Delta = min(scan_range, window - word_length)` diagonals
   (`-off_diagonal_range`, default 0, :657, :686-716).
4. Ungapped extension (:723-740): exact matrix extension
   `s_NuclUngappedExtendExact` (:149-253) when `matrix_only_scoring` (rmblastn)
   or `word_length < 11`; otherwise the fast 4-base-table version
   `s_NuclUngappedExtend` (:262-357) which re-runs exactly if the approximate
   score >= `reduced_nucl_cutoff_score` (0.8 x cutoff, `blast_parameters.c:410-414`).
   X-drop: `-cutoffs->x_dropoff` (raw).
5. Saved if `score >= cutoffs->cutoff_score` (or an off-diagonal partner was
   found) via `BLAST_SaveInitialHit` (:742-752); the diagonal's `last_hit` is
   advanced to the end of the ungapped alignment (:754-756).

## 2. Ungapped cutoff and X-drop (`blast_parameters.c`)

- `x_dropoff_init = ceil(x_dropoff_bits * ln2 / Lambda_ungapped)` per context
  (:213-215). `-xdrop_ungap` default 20 bits (`BLAST_UNGAPPED_X_DROPOFF_NUCL`).
- `BlastInitialWordParametersUpdate` (:280-419): `gap_trigger` raw =
  `(27 bits * ln2 + logK) / Lambda` with **ungapped** K-A values (:336-346,
  `BLAST_GAP_TRIGGER_NUCL` 27.0, API `SetGapTrigger`, no blastn CLI flag).
  Gapped search: `cutoff_score = min(gap_trigger, hit_params cutoff_score_max)`
  (:367-373). Ungapped search (`-ungapped`): `BLAST_Cutoffs` with
  `CUTOFF_E_BLASTN` 0.05 (`blast_parameters.h:76`) on search space
  `min(2*qlen, subj) * subj` (:348-366); blastn does **not** cap at gap_trigger (:364-365).
- Notes (arithmetic, not executed): blastn 2/-3 ungapped Lambda 0.55, K 0.21
  (`blast_stat.c:674-675`) -> gap_trigger ~31 raw, X-drop ~26 raw;
  megablast 1/-2 (Lambda ~1.28, K ~0.46, :660-661) -> ~14 raw, X-drop 11 raw;
  blastn-short 1/-3 (1.374, 0.711, :638-639) -> ~13 raw.

## 3. Preliminary gapped extension: `BLAST_GetGappedScore` (`blast_gapalign.c:3664-4121`)

1. Ungapped HSPs arrive sorted by score. Each is skipped if an interval tree
   already holds a higher-scoring HSP enveloping it on both ends within
   `min_diag_separation` diagonals (:3907-3908; `blast_itree.c:810-847`;
   50 for blastn, 6 for megablast, `blast_nucl_options.cpp:239,259`).
2. `cutoff = hit_params->cutoffs[context].cutoff_score` (:3917-3920) =
   raw score for `-evalue` at the context's effective search space.
3. Greedy (`ePrelimGapExt == eGreedyScoreOnly`, megablast default): seed at
   the middle of the ungapped HSP (:4019-4022), `BLAST_GreedyGappedAlignment`
   (:2762-2947) with `X = gap_align->gap_x_dropoff` (raw of 25 bits,
   `BLAST_GAP_X_DROPOFF_GREEDY`), traceback kept. Greedy memory sized by
   `max_dist = min(GREEDY_MAX_COST 1000, subject_len/2 + 1)` (:322-331,
   `greedy_align.h:47-50`); non-affine 0/0 costs use `gap_extend = reward/2 + penalty`
   (:199-200); odd rewards are doubled internally (:183-189).
4. DP (`eDynProgScoreOnly`, blastn/dc-megablast or `-no_greedy`,
   `blast_args.cpp:664-667,684-688`): `s_BlastDynProgNtGappedAlignment`
   (:2949-3032) -> `s_BlastAlignPackedNucl` (:3034) left then right from a
   4-base-aligned start; X-drop is capped at the ungapped score (:2963-2965);
   `-xdrop_gap` default 30 bits (`BLAST_GAP_X_DROPOFF_NUCL`).
5. `gap_align->score >= cutoff` -> `Blast_HSPInit` + `Blast_HSPListSaveHSP`
   + interval tree (:4068-4100), else dropped (greedy traceback freed :4103-4107).
6. Raw X-drop conversion: `gap_x_dropoff = bits*ln2/min_gapped_lambda`,
   `gap_x_dropoff_final = max(final_bits*ln2/lambda, gap_x_dropoff)`
   (`blast_parameters.c:455-461`). `-xdrop_gap_final` default 100 bits
   (`BLAST_GAP_X_DROPOFF_FINAL_NUCL`), applied in traceback
   (`blast_traceback.c:1470, 1424-1433`).
7. Ungapped-only run (`-ungapped`): `BLAST_GetUngappedHSPList` (:4723) and
   later `Blast_HSPListReevaluateUngapped` + linking (`blast_engine.c:1489-1544`).

## 4. Traceback stage differences (`blast_traceback.c:259-731`)

Greedy traceback for megablast (:503-506) vs `ALIGN_EX` affine traceback
(:507-512, `blast_gapalign.c:4553-4719`, prunes leading/trailing gaps
:4670-4705). Start point from `BlastGetStartForGappedAlignmentNucl`
(`blast_gapalign.c:3323-3391`). Afterwards HSPs are rescored with ambiguities
(`Blast_HSPReevaluateWithAmbiguitiesGapped`, `blast_hits.c:479`, cutoff =
`hit_params->cutoffs[context].cutoff_score` :513), filtered by
identity/length, and purged for containment.

## 5. HSP saving and per-subject limits (`blast_hits.c`)

- `Blast_HSPListSaveHSP` (:1754-1809): array grows to `hsp_max`
  (`BlastHspNumMax` :213-229: `hsp_num_max` or INT4_MAX by default); when
  full it becomes a min-heap by score and worse HSPs are discarded.
- Chunk merge: `Blast_HSPListsMerge` (:2857) with `hsp_num_max`.
- `Blast_HSPListPurgeHSPsWithCommonEndpoints` (:2455-2535): HSPs sharing a
  start or end point (same context/frame) are deleted (prelim :549 in
  `blast_engine.c`; traceback :650, :682).
- `Blast_HSPListGetEvalues` (:1811-1905), reap by `prelim_evalue`
  (`blast_engine.c:643-675`) and finally by `-evalue` (`Blast_HSPListReapByEvalue` :1976).

## 6. Hit list size, culling, best hit, subject best hit

- `hitlist_size` = number of subject sequences kept per query
  (`-max_target_seqs`, default 100 `cmdline_flags.cpp:91`; set in
  `blast_args.cpp:2921-2975`; API default 500). Prelim keeps
  `GetPrelimHitlistSize` = `min(max(2h,10), h+50)` for gapped (`blast_hits.c:44-72`).
  `Blast_HitListUpdate` (:3243-3297) heapifies by best E-value once full and
  drops the worst subject. `s_BlastPruneExtraHits` (`blast_traceback.c:878-893`)
  trims to `hitlist_size` after traceback. Warning if < 5 (`blast_args.cpp:2973`).
  For a genome in few chromosomes this limit is on *chromosomes*, not HSPs.
- `-max_hsps N` -> `max_hsps_per_subject`; applied post-traceback by
  `Blast_TrimHSPListByMaxHsps` (`blast_hits.c:2049-2069`) keeping the N best
  by score (`s_FilterBlastResults` `blast_traceback.c:837-867`).
- `-culling_limit N` (`blast_args.cpp:3291-3299`) -> culling writer
  (`setup_factory.cpp:333-342`, prelim limit N+3). `hspfilter_culling.c`:
  an HSP is deleted once N higher "dominating" HSPs overlap > 50% of its
  query range (`s_DominateTest` :79-120: `d = 4*s1*l1 + 2*s1*l2 - 2*s2*l1 - 4*s2*l2 >= 0`).
- `-best_hit_overhang` / `-best_hit_score_edge` (defaults 0.1/0.1,
  `hspfilter_besthit.h:86-98`; algorithm described :50-82; prelim overhang
  forced to 0 `setup_factory.cpp:327`).
- `-subject_besthit` (`blast_args.cpp:3327`): `Blast_HSPListSubjectBestHit`
  (`blast_hits.c:2538-2606`) removes, within one subject, lower-scoring HSPs
  whose query range lies within a higher-scoring HSP's range +/- 3
  (`DEFAULT_SUBJECT_BESTHIT_NUCL_MAX_RANGE_DIFF`, `blast_options.h:1212`),
  including the opposite strand for blastn (:2581-2602). Applied per chunk in
  prelim (`blast_engine.c:586-591`) and after traceback (`blast_traceback.c:858-862`).
  For a repeat consensus vs a chromosome this keeps only the best copy per
  query region: **do not use for repeat annotation**.

## 7. Identity / coverage / length filters

- `s_HSPTest` (`blast_hits.c:993-1001`): delete if
  `num_ident*100 < align_length*percent_identity` or `align_length < min_hit_length`.
  Called only during traceback (`blast_traceback.c:606-620, 676-678`), so
  `-perc_identity` never reduces prelim work.
- `Blast_HSPGetQueryCoverage` (:1034-1043) = `100*(qend-qstart)/qlen` (+0.5 if
  < 99); `-qcov_hsp_perc` applied post-traceback (`blast_traceback.c:1764`,
  :837-867) and in the ungapped path (`blast_engine.c:1538`).
- No minimum HSP length other than `min_hit_length` (no CLI) and the word size.

## 8. Sum statistics / `link_hsps.c`

`do_sum_stats` is FALSE for gapped blastn, TRUE for `-ungapped`
(`blast_options.c:1464-1470`); `-sum_stats` is hidden for blastn
(`blast_args.cpp:128,240-245`). When on, `BLAST_LinkHsps`
(`link_hsps.c:1761-1811`) does even-gap linking (`s_BlastEvenGapLinkHSPs`
:415) with `gap_prob 0.5`, `gap_decay_rate 0.5`, `gap_size 40`,
`overlap_size 9` (`blast_parameters.h:66-71`, `BlastLinkHSPParametersNew`
`blast_parameters.c:598-622`); `cutoff_small_gap` = ungapped cutoff
(:625-648). So for default gapped blastn each HSP gets its own E-value.
