# Control flow

All line numbers at commit `3c28777`.

## A. Program start and option handling (`main.c`)

1. `mm_set_opt(0, ...)` loads defaults (`main.c:148` -> `options.c:93-95`).
2. First pass over argv applies only `-x PRESET` (`main.c:150-163` ->
   `mm_set_opt(preset)` `options.c:91-193`). Unknown preset -> exit 1.
3. Second pass applies every other option in command-line order
   (`main.c:166-351`). Notable conversions: `-f` float vs int (`main.c:325-331`),
   `-U lo,hi` (`main.c:322-324`), `-r bw[,bw_long]` (`main.c:319-321`),
   `-z zdrop[,zdrop_inv]` (`main.c:341-343`), `-O/-E` dual values
   (`main.c:344-349`), `--secondary=yes|no` toggles `MM_F_NO_PRINT_2ND`
   (`main.c:282-283`), `--rmq[=yes|no]` (`main.c:301-303`), `--sr[=dna|rna|no]`
   (`main.c:269-279`), `-N` remembers the previous `best_n` (`main.c:176`).
4. If neither `-d` nor CIGAR output (`-c`/`-a`) is requested, the index is
   built without sequence (`MM_I_NO_SEQ`, `main.c:352-353`).
5. `mm_check_opt` validates (`main.c:354`, `options.c:202-277`): `-k/-w > 0`,
   `0 <= -p <= 1`, dual gap constraint `E1 > E2 && O1+E1 < O2+E2`
   (`options.c:246-250`), `(O1+E1)+(O2+E2) <= 127` (int8 DP), `--score-N < -B`,
   `zdrop >= zdrop_inv`, `-P` incompatible with `--secondary=no`, `--rmq`
   incompatible with `--sr`/`--splice`.
6. `-N 0` is rewritten to the previous `-N` plus `--secondary=no` with a
   warning (`main.c:356-359`).
7. `mm_idx_reader_open` (`main.c:425`, `index.c:595-611`) detects `.mmi` by
   magic `MMI\2` (`index.c:571-593`).
8. Loop `while ((mi = mm_idx_reader_read(...)))` (`main.c:437-506`): one
   iteration per index part (`-I`). Per part: SAM header (`main.c:445-461`),
   `mm_mapopt_update` derives `mid_occ` (`main.c:465`, see C.0), optional
   BED/junction/alt loading, then `mm_map_file` or `mm_map_file_frag` over all
   query files (`main.c:493-500`). The query files are re-read for every part.
9. `--split-prefix` merge (`main.c:510-511`), flush, timing/RSS summary
   (`main.c:518-524`).

## B. Index construction (`index.c`)

1. `mm_idx_gen` (`index.c:389-408`) creates `mm_idx_t` with `1<<b` empty
   buckets (`mm_idx_init` `index.c:50-60`) and runs a 3-step `kt_pipeline` on
   `min(n_threads, 3)` threads (`index.c:399`).
2. Step 0 reads a mini-batch of `mini_batch_size` (50 Mbp) reference bases
   (`index.c:323`), stops when `sum_len > batch_size` (`-I`, `index.c:321`),
   grows `seq[]`, copies names and 4-bit sequence into `S` unless `MM_I_NO_SEQ`
   (`index.c:344-366`).
3. Step 1 sketches each sequence: `mm_sketch(..., w, k, rid, is_hpc, &s->a)`
   (`index.c:374`).
4. Step 2 appends each minimizer to the bucket chosen by the low `b` bits of
   its hash (`mm_idx_add`, `index.c:306-313`; 16 B each).
5. `mm_idx_post` -> `kt_for(worker_post)` over buckets (`index.c:280-283`,
   `226-278`): radix sort the bucket by minimizer, count distinct keys, build
   khash sized to `n_keys`, copy positions of multi-occurrence minimizers into
   `p[]` sorted by position, free the temp array.
6. If `-d`, `mm_idx_dump` writes the part (`index.c:631`, format in
   `io-formats.md`).

### B.1 `mm_sketch` (`sketch.c:77-143`)

- Walks the sequence keeping forward and reverse-complement k-mers
  (`sketch.c:106-107`); ambiguous bases reset the k-mer (`sketch.c:115`).
- HPC (`-H`): each homopolymer run counts as one base; `kmer_span` is the
  actual span in bp of the k HPC-bases (`sketch.c:94-104`); non-HPC span = k.
- Skips k-mers equal to their reverse complement (`sketch.c:108`); strand bit
  = which of the two is smaller (`sketch.c:109`).
- `info.x = hash64(canonical_kmer) << 8 | span` (`sketch.c:112`). The window
  of `w` consecutive k-mers emits its minimum; ties (identical k-mers) are all
  emitted (`sketch.c:117-137`). With `w = 1` every k-mer is emitted.

## C. One query through `mm_map_frag_core` (`map.c:227-378`)

0. (per index part, before any query) `mm_mapopt_update` (`options.c:69-83`):
   if `mid_occ <= 0`, `mid_occ = mm_idx_cal_max_occ(mi, mid_occ_frac)` = the
   occurrence count at the `(1 - f)` quantile over *distinct* minimizers, +1
   (`index.c:198-220`; `f <= 0` -> `INT32_MAX`), then clamped to
   `[min_mid_occ, max_mid_occ]` (`options.c:75-78`). With `-v3` it prints
   `mid_occ = N` (`options.c:81-82`).
1. Guard: empty / too many segments / `--max-qlen` (`map.c:243-244`). Per-query
   hash from name, length and `--seed` for deterministic tie-breaking
   (`map.c:246-248`).
2. `collect_minimizers` sketches the query with the index's `w,k,HPC`
   (`map.c:59-72`, `mm_sketch` call `map.c:65`); `-T` SDUST masking if > 0
   (`map.c:68-69`, default off).
3. `mm_seed_mz_flt` (`map.c:251`, `seed.c:5-28`): drop a query minimizer that
   occurs more than `mid_occ` times *within the query* and more than
   `q_occ_frac` (1%) of all query minimizers. Irrelevant for short queries.
4. `collect_seed_hits` (default) or `collect_seed_hits_heap` (`MM_F_HEAP_SORT`,
   set by `sr`) (`map.c:252-253`) -> `mm_collect_matches` (`seed.c:98-132`):
   a. `mm_seed_collect_all` looks each query minimizer up with `mm_idx_get`
      (`seed.c:42`); misses are skipped; `is_tandem` marked (`seed.c:47-48`).
   b. If `occ_dist > 0 && max_max_occ > max_occ`: `mm_seed_select`
      (`seed.c:56-96`) marks every seed with `n > max_occ` as filtered, then
      for each maximal run of such seeds un-filters the `round(run_bp /
      occ_dist)` lowest-occurrence ones (<= 128, `seed.c:54,73-76`), re-filtering
      any with `n > max_max_occ` (`seed.c:89-91`). Otherwise all seeds with
      `n > max_occ` are filtered (`seed.c:109-111`). Here `max_occ` is
      `opt->mid_occ` on the first pass (`map.c:252-253`).
   c. Filtered seeds contribute their query span to `rep_len` (`seed.c:117-122`);
      kept seeds go to `mini_pos[]` and `n_a += n` (`seed.c:124-126`).
      `--dbg-seed-occ` prints `SF pos occ flt` per seed (`seed.c:115-116`).
   d. Anchors are materialised for every occurrence of every kept seed
      (`map.c:176-200`), `skip_seed` applies `-D/-X/--for-only/--rev-only`
      (`map.c:78-100`), then radix sort (`map.c:202`). `--print-seeds` dumps
      `SD` lines (`map.c:255-260`).
5. Chaining limits (`map.c:262-274`): query gap = `max_gap` (or `max(qlen,
   max_gap)` for `--sr`), ref gap = `max_gap_ref` if set, else derived from
   `max_frag_len`, else `max_gap`; `chn_pen_gap = chain_gap_scale * 0.01 * k`,
   `chn_pen_skip = chain_skip_scale * 0.01 * k` (0 by default).
6. Chain: `mg_lchain_rmq` if `--rmq` else `mg_lchain_dp` (`map.c:275-281`,
   see C.1). Output: `n_regs0` chains, `u[]`, compacted `a[]`.
7. Re-chaining (`map.c:283-316`): (a) long-read "long join": if `bw_long > bw`
   and the best chain leaves >= `rmq_rescue_size` (1000 bp) unaligned or
   covers > `rmq_rescue_ratio` of the query, re-chain all chained anchors with
   `mg_lchain_rmq` at `bw_long` (`map.c:283-292`); (b) `--sr` path: if
   `max_occ > mid_occ` and some seeds were filtered and the best chain does not
   cover all segments (or there is none), redo seeding with `max_occ` instead
   of `mid_occ` and re-chain (`map.c:293-316`). Only `sr` sets `max_occ`
   (5000, `options.c:167`); `-f INT1,INT2` also sets it (`main.c:331`).
8. `mm_gen_regs` converts chains to `mm_reg1_t` sorted by score desc with hash
   tie-break (`map.c:320`, `hit.c:52-88`). ALT-contig handling if `--alt`.
9. `chain_post` (`map.c:332`, `206-213`): unless `-P`, `mm_set_parent` marks
   secondaries by query overlap (`hit.c:125-186`; test at `hit.c:165`: overlap
   fraction of the shorter chain minus uncovered fraction > `mask_level` and
   uncovered <= `mask_len`), then `mm_select_sub` keeps primaries plus at most
   `best_n` secondaries with `score >= pri_ratio * parent.score || score + 2k
   >= parent.score` (`hit.c:265-267`), plus opposite-strand secondaries with
   score > `0.8 * max_gap` (`hit.c:268-271`, `map.c:210`). `pri_ratio == 0`
   disables the whole function (`hit.c:257`).
10. `mm_est_err` computes `dv` from the fraction of query minimizers that are
    anchors in the chain (`map.c:334`, `esterr.c:30-64`), and strand-retained
    secondaries with much higher divergence are dropped (`map.c:335`,
    `hit.c:283-299`). Skipped for `--sr`.
11. `align_regs` (`map.c:339`, `215-225`): returns immediately without CIGAR
    (`map.c:217`); else `mm_align_skeleton` (C.2) then `mm_set_parent` /
    `mm_select_sub` (check_strand = 0) / `mm_set_sam_pri` again
    (`map.c:219-223`).
12. `mm_set_mapq2` (`map.c:341`, `hit.c:432-485`): MAPQ from chain score,
    `subsc`, `dp_max/dp_max2`, identity, anchor count, and `uniq_ratio =
    sum_sc / (sum_sc + rep_len)` (`hit.c:445`) so filtered repeat seeds lower
    MAPQ; minus `4.343*ln(n_sub+1)`; capped at 60.
13. Multi-segment (paired) queries split chains per segment (`mm_seg_gen`) and
    pair (`mm_pair`) (`map.c:343-355`). Arena reset check (`map.c:366-377`).

### C.1 `mg_lchain_dp` (`lchain.c:148-217`)

- Arrays `p` (predecessor, 8 B), `f` (score, 4 B), `v` (peak, 4 B), `t`
  (visit marks, 4 B) of size `n` anchors (`lchain.c:163-166`): 20 B/anchor.
- `max_dist_x = max(max_dist_x, bw)`, same for `y` (`lchain.c:160-161`);
  `max_drop = bw` (`lchain.c:151`).
- For each anchor `i`: `st` advances so that `a[st]` is on the same strand/rid
  and within `max_dist_x` on the reference (`lchain.c:172`); at most `max_iter`
  (5000) predecessors are scanned (`lchain.c:173`).
- `comput_sc(ai, aj)` (`lchain.c:113-138`): `dq` = query gap, `dr` = ref gap;
  reject if `dq <= 0`, `dq > max_dist_x`, `dr == 0`, `dq > max_dist_y`, or
  `|dr - dq| > bw`. Score = `min(q_span_j, min(dr, dq))` minus, when there is
  a diagonal shift or the gap exceeds the span, `chn_pen_gap*|dd| +
  chn_pen_skip*dg + 0.5*log2(|dd|+1)` (`lchain.c:124-135`).
- Early exit when `max_skip` (25) predecessors already on a better chain are
  seen (`lchain.c:182-184`); a "best in window" shortcut `max_ii`
  (`lchain.c:189-204`) keeps long chains linked beyond `max_iter`.
- `mg_chain_backtrack` (`lchain.c:27-76`): candidate chain ends are anchors
  with `f >= min_sc`, sorted by score; each is traced back through unused
  anchors until the score drop exceeds `max_drop`; kept only if `score >=
  min_sc && n_anchors >= min_cnt`. Chains are then compacted (`lchain.c:216`).

### C.2 `mm_align_skeleton` / `mm_align1` (`align.c:1048-1120`, `645-914`)

1. Encode query both strands (`align.c:1056-1061`); drop anchors of discarded
   chains (`mm_squeeze_a`, `align.c:1064`).
2. For each hit `mm_align1` (`align.c:1098`); splice modes align twice
   (`align.c:1068-1096`).
3. `mm_align1`: build scoring matrix from `-A/-B/-b(transition)/--score-N`
   (`align.c:659`, `ksw_gen_ts_mat` `align.c:26-36`); alignment bandwidth =
   `1.5 * bw + 1` (`align.c:660-661`).
   - Non-`sr`: `mm_fix_bad_ends` trims chain ends with large indels or until
     `2 * min_chain_score` matched bases (`align.c:675`, `527-561`);
     `mm_filter_bad_seeds` marks seeds implying both a long insertion and a
     long deletion as `MM_SEED_IGNORE` (`align.c:677`, `454-489`);
     `mm_filter_bad_seeds_alt` marks long-join seeds (`align.c:678`, `491-525`).
   - DP window `[rs0,re0] x [qs0,qe0]`: from first/last anchor, extended
     outward by the unaligned query length (capped at `max_gap`) and by
     `(l*a - q)/e` extra reference bases for possible gaps (`align.c:704-759`).
     `sr`: whole query, reference padded by `(l*a + end_bonus - q)/e`
     (`align.c:696-703`).
4. Left extension with `ksw_extz2_sse`/`ksw_extd2_sse` in `EXTZ_ONLY|RIGHT|
   REV_CIGAR` mode, `end_bonus`, `zdrop` (`align.c:779-799`, `mm_align_pair`
   `336-368`). `mm_align_pair` skips DP (marks `zdropped`) when `tlen*qlen >
   max_sw_mat` (1e8 cells, `align.c:349-351`). Dual-affine `extd2` is used
   when `q != q2 || e != e2` (`align.c:356-359`), i.e. by default.
5. Gap filling between consecutive kept anchors, but a DP is launched only at
   the last anchor, at a `LONG_JOIN` anchor, or when both the query and
   reference gaps since the last DP are >= `min_ksw_len` (200) (`align.c:810`).
   `sr`: ungapped scoring first, DP only if it scores below the best possible
   gapped score (`align.c:823-834`). Normal: first pass with `KSW_EZ_APPROX_MAX`
   (`align.c:840`), then `mm_test_zdrop` on the CIGAR; if a drop > `zdrop`
   (or > `zdrop_inv` with inversion signature) is found, re-align exactly
   (`align.c:843-844`). If `ez->zdropped`, the chain is split at the drop into
   `r2` when the remainder has >= `min_cnt` anchors (`align.c:848-868`), and
   `mm_align1_inv` later tries an inverted alignment between the two pieces
   (`align.c:1103-1108`).
6. Right extension (`align.c:874-890`); final coordinates (`align.c:893-895`);
   `mm_update_extra` recomputes `mlen/blen/dp_max/n_ambi` from the CIGAR
   (`align.c:906`, `254`).
7. Back in the skeleton: `mm_filter_regs` (`align.c:1113`, `hit.c:301-320`)
   drops `cnt < min_cnt`, `mlen < min_chain_score`, `dp_max < min_dp_max`, or
   both ends clipped > `max_clip_ratio` (1.0 = never). For queries >=
   `rank_min_len` (500 bp) and not `sr`/`-P`, `mm_update_dp_max` rescales
   `dp_max` of all hits by a divergence-aware mismatch cost when the second
   best is within `rank_frac` (0.9) of the best (`align.c:1114-1117`,
   `1022-1046`). `mm_hit_sort` orders by `dp_max` (`align.c:1118`, `hit.c:188-218`).

## D. Output step (`map.c:578-641`)

For each query, each kept `mm_reg1_t` is written unless `MM_F_NO_PRINT_2ND`
and it is secondary (`map.c:612-613`); PAF via `mm_write_paf4` (`map.c:617`),
SAM via `mm_write_sam3` (`map.c:615`). Unmapped records only with
`--paf-no-hit` or SAM without `--sam-hit-only` (`map.c:620-626`). Everything is
freed per batch (`map.c:628-636`).

## Threading / parallelism model

- Mapping: `kt_pipeline(pl_threads, worker_pipeline, 3 steps)` with
  `pl_threads = 1 if -t 1 else (3 if -2 else 2)` (`map.c:677-678`,
  `kthread.c:130-159`). Step 0 reads `-K` bases (`bseq.c:80-119`), step 1 maps
  with `kt_for(n_threads, worker_for, n_frag)` (`map.c:576`), step 2 writes in
  input order. `kt_for` hands out indices `i, i+n, i+2n, ...` per thread and
  then work-steals from the most-behind thread (`kthread.c:30-52`), so one
  slow query (huge repeat) delays only its batch's completion.
- Indexing: 3-step pipeline on <= 3 threads (`index.c:399`), then `kt_for`
  over `1<<b` buckets for sorting/hashing (`index.c:282`).
- Each mapping thread has its own `mm_tbuf_t` arena; the index is read-only
  and shared. Memory per thread while mapping a query ~ `36 B * n_anchors`
  (16 B anchors + 20 B chaining arrays, `lchain.c:163-166`) plus DP buffers.

## Error handling conventions

- Fatal CLI/config errors print `[ERROR]` and `return 1` (`main.c:153-161,
  354-355, 426-434`). Warnings `[WARNING]` at `mm_verbose >= 2`; progress
  `[M::func::time*cpu]` at `>= 3` (default 3, `main.c:145`).
- Library functions return `0`/`NULL` on failure (`mm_idx_reader_read`,
  `mm_map_file_frag` returns -1 on unopenable input, `map.c:667-671`).
- `assert()` is used liberally for invariants (compiled in by default,
  `CFLAGS` has no `-DNDEBUG`, `Makefile:1`), e.g. HPC with `--sr`
  (`align.c:655`), arena leak check (`map.c:370`).
