# Control flow

## `ssearch36 -n [opts] query.fa genome.fa` (threaded build)

1. `main()` (`src/comp_lib9.c`) calls `initenv()` (`src/doinit.c:~330-500`):
   generic options (`-m -N -T -Z -w -O -q ...`) are handled in `doinit.c:376-490`,
   algorithm options are forwarded to `f_getopt()` (`src/initfa.c:1030-1327`) and
   `-X...` to `parse_ext_opts()` (`initfa.c:1346-1440`). `-n` sets
   `qdnaseq=SEQT_DNA` (`initfa.c:1227-1232`).
2. Query is read with `QGETLIB` into `aa0[0]`, max `MAXTST` 40,000 residues
   (`comp_lib9.c:655-656,710-715`; `defs.h:92`). Longer queries are truncated with
   a warning.
3. `resetp()` (`initfa.c:1570-1747`) switches to DNA: alphabet, `e_cut/=5`
   (10 -> 2.0) and `e_cut_r = 0.001` (`1584-1587`), both strands (`qframe=2`,
   `1583`), gaps -12/-4 (`1680-1686`), matrix `+5/-4` built by `mk_n_pam()`
   (`1729-1737`). `last_init()` computes ungapped lambda/K (`initfa.c:1823+`).
4. Chunk geometry: `max_tot = MAXLIB` (150,000, `initfa.c:603`, `defs.h:93`);
   `maxn = reset_maxn(&m_msg, 150, max_tot)` (`comp_lib9.c:894`,
   `compacc2.c:155-174`): `-N` lowers `maxn`; `l_overlap = 150`;
   `maxt3 = maxn - 150`.
5. Buffer sizing (`comp_lib9.c:805-848`): `max_chain_seqs = 6,000,000/1200 = 5000`
   entries per buffer for DNA (min 600); `max_work_buf` between 2x and 4x
   `fa_max_workers`; each buffer's sequence area is
   `max(max_chain_seqs*ave_seq_len, 2*max_tot)` = 6 MB (`comp_lib9.c:2039`).
   `alloc_comp_bufs()` allocates them (`comp_lib9.c:1983+`).
6. `init_thr()` (`src/pthr_subs2.c:48-120`) creates `fa_max_workers` pthreads
   running `work_thread()` (`src/work_thr2.c:108`). Each worker copies `pstruct`
   and `pam2` (`work_thr2.c:255-279`), builds its own forward and
   reverse-complement query profiles via `init_aa0()` -> `init_work()`
   (`work_thr2.c:313-315`; `dropgsw2.c:132-620`).
7. Reader loop `seqr_chain_work()` -> `next_seqr_chain()` (`comp_lib9.c:2500-2662`):
   `GETLIB` (= `agetlib()` for FASTA files, `src/nmgetlib.c:574-705`) reads up to
   `maxt` residues of the current entry; if the entry continues, `lcont` is set
   (`nmgetlib.c:699-701`), the last 150 residues are saved (`comp_lib9.c:2644-2646`)
   and prepended to the next chunk (`2505-2508`), `cont` increments (`2546`),
   `l_offset += n1 - 150` (`2655`). Continuation chunks are not counted as new
   `db.entries` (`2616-2620`). Filled buffers go to workers via
   `put_rbuf()`/`get_wbuf()` (`pthr_subs2.c:188,322`).
8. Worker: `buf_do_work()` (`src/compacc2.c:2960-3057`) loops over entries,
   skips those outside `-M` (`2990-2996`), calls `do_work()` for each
   (chunk, strand) (`3012-3015`). `do_work()` (`dropgsw2.c:696-870`) runs
   `smith_waterman_sse2_byte()`, re-runs `_word()` on overflow (`801-826`), and
   stores the single best score in `rst->score[0]` (`829`). With `-z 11..16`
   `buf_shuf_work()` also scores a shuffled copy (`work_thr2.c:357-359`).
9. Main thread: `save_best()` (`compacc2.c:2369-2616`) per returned buffer:
   accumulates `ldb` counts, appends the best-over-strands score of each chunk to
   `stats[]` (<= 60,000; `2496-2545`), and keeps every (chunk, strand) with
   `zscore > zbestcut` as a `beststr` (`2558-2611`; `zbestcut` starts at
   `-FLT_MAX`, `comp_lib9.c:917`, so until 60,000 hits everything is kept).
   When `stats[]` fills, `process_hist()` is run early and later scores are
   z-scored on the fly (`2500-2516`).
10. After the scan: `process_hist()` fits statistics (`comp_lib9.c:1048-1070`;
    `scaleswn.c:236-369`); if `nbest < shuff_max` (500) and `-z` is not 3/-1,
    all hits are shuffled with `buf_shuf_seq()` and refitted as `-z 1x`
    (`comp_lib9.c:1079-1096`); `zdb_size` defaults to `ldb.entries` (`1103`);
    `-z 2x` shuffles the top hits (`1174-1216`). `scale_scores()` assigns
    E-values and sorts (`comp_lib9.c:1165`; `scaleswn.c:2942-2967`).
11. `nshow` = number of hits with `escore < e_cut` after skipping `-F e_low`
    (`comp_lib9.c:1340-1391`), modified by `-b`. `-E` and `-z -1` together mean
    `nshow = min(mshow, nbest)` with `mshow` default 50 (`1384-1385`,
    `doinit.c:296`).
12. Alignment pre-pass (needed for `-m 8/9/B`): `repeat_thresh` per hit
    (`comp_lib9.c:1533-1550`), then `buf_align_seq()` farms `do_walign()` to the
    workers (`1554`; `compacc2.c:3060` `buf_do_align()`). `do_walign()`
    (`dropgsw2.c:1001-1040`) calls `nsw_malign()` (`src/wm_align.c:149-297`) which
    aligns with scalar `sw_walign()` (`wm_align.c:37`), then recurses on the
    library prefix and suffix while scores exceed `repeat_thresh`.
13. Output: `showbest()` (`src/mshowbest.c:112-770`) prints one line per
    alignment in each hit's `a_res` chain (`do { ... } while (cur_ares_p = cur_ares_p->next)`,
    `504-741`); `showalign()` (`src/mshowalign2.c:137+`) prints each alignment
    block (`while (cur_ares_p != NULL && cur_ares_p->nres > 0)`, `485`).
14. Multi-query: the loop returns to step 2 (`comp_lib9.c:1830-1900`); the
    library stays in memory if it fit in `max_memK` (`1889-1899`).

## Threading / parallelism model

- Producer/consumer with a fixed ring of buffers: main thread reads and fills,
  `fa_max_workers` threads score; mutex+condvar pairs `reader_mutex`/
  `worker_mutex` (`pthr_subs2.c:76-93`). A buffer is the unit of work, not a
  sequence; `buf_do_work()` processes all entries in it serially.
- `-T n` sets `fa_max_workers` (`doinit.c:432-444`); default = number of CPUs
  (`sysconf(_SC_NPROCESSORS_CONF)`, `doinit.c:255`). `-DMAX_WORKERS=8` in the
  Makefile only initialises the variable (`doinit.c:92`).
- `-DCOMP_THR` selects the threaded driver (`make/Makefile.pcom_t:128-129`);
  lalign36 links the non-threaded `comp_mlib9.o` (`Makefile.pcom_t:50`).
- Alignment phase is also distributed to workers (`buf_align_seq()`), shuffles
  too (`BUF2_DOSHUF`).
- No worker-side filtering: every chunk's score goes back to the main thread.

## Error handling conventions

`fprintf(stderr, "*** ERROR [%s:%d] ...", __FILE__, __LINE__)` then `exit(1)`
for allocation and size failures; `*** Warning` for recoverable issues (e.g.
query truncation `comp_lib9.c:712`). Matrix values that do not fit 8-bit lanes
are fatal (`dropgsw2.c:596-601`).
