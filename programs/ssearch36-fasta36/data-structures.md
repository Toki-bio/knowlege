---
title: "ssearch36-fasta36: Core data structures"
type: reference
generated: true
---

> Part of the [ssearch36-fasta36](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/ssearch36-fasta36/references/data-structures.md)

{% raw %}
# Core data structures

## `struct pstruct` (`src/param.h:50-131`) -- all scoring/statistics parameters

| Field | Meaning | Set by |
|---|---|---|
| `gdelval`, `ggapval` | gap open / extend (negative) | `-f`, `-g` (`initfa.c:1143-1151`); DNA defaults -12/-4 (`initfa.c:1680-1686`) |
| `p_d_mat`, `p_d_mis`, `p_d_set` | DNA match/mismatch | `-r "+5/-4"` (`initfa.c:1248-1258`); default 5/-4 (`638-639`) |
| `pam2[2]`, `pam2p[2]` | scoring matrix [0]=plain, [1]=lower-case-as-X (`-S`) | `mk_n_pam()` / `standard_pam()` / `initpam()` |
| `pam_l`, `pam_h` | lowest / highest matrix value; `pam_l` drives x-drop in fasta | |
| `zsflag`, `zsflag2` | `-z` method and the method for `-z 2x` shuffles | `initfa.c:1312-1325` |
| `e_cut`, `e_cut_r`, `do_rep` | display E() cutoff; cutoff for additional alignments; enable sub-alignments | `-E` (`initfa.c:1113-1141`) |
| `zdb_size`, `zdb_size_set` | `-Z` database size | `doinit.c:488-490` |
| `n1_low`, `n1_high` | `-M` length window (applies to chunk length) | `initfa.c:1211-1226` |
| `sw_flag` | 1 = Smith-Waterman alignments; 0 = banded (fasta DNA) | `-A` / `-XA` (`initfa.c:1055-1058,1365-1369`) |
| `max_repeat`, `repeat_thresh` | lalign `-K` (default 50 `doinit.c:288`) | |
| `pLambda, pK, pH`, `LK_set` | ungapped Karlin-Altschul params used for fasta thresholds and `-z 3` | `last_init()` |
| `param_u.fa` (`struct fastr`, `param.h:14-36`) | `ktup`, `optcut`, `optwid`, `cgap`, `E_join`, `E_band_opt`, `use_E_thresholds`, `optflag` | fasta36 only |

Shared read-only by the main thread; each worker copies it (`work_thr2.c:255-279`).

## `struct mngmsg` (`src/structs.h:75-191`) -- run-wide "manager message"

Key fields: `n0` query length; `qdnaseq`; `qframe` (2 = both strands for DNA,
`initfa.c:1583`); `nframe`/`nitt1`; `e_cut`, `e_low` (`-F`); `mshow`/`mshow_min`
(`-b`), `ashow` (`-d`), `nshow` (computed count below `-E`, `comp_lib9.c:1364-1385`);
`markx`, `show_code`, `markx_list` (`-m`); `shuff_max` (`-k`, default 500
`initfa.c:969`); `max_memK` (`-XM`); `ldb_info` (`struct lib_seq_info`,
`structs.h:33-40`: `maxn` chunk size, `maxt3` continuation size, `l_overlap`=150,
`ldnaseq`, `term_code`); `db` (residues/entries read) and `ldb` (counted in
save_best); `pstat_void`/`pstat_void2` (fitted statistics); `aln` (display widths).

## `struct seq_record` / `struct mseq_record` (`src/param.h:176-206`) -- one library chunk

`aa1b` encoded residues (NUL-terminated, padded `SEQ_PAD` 16); `n1` chunk length
(includes the 150 overlap for continuations, `comp_lib9.c:2619`); `l_offset`
0-based residues skipped before this chunk (`comp_lib9.c:2547`); `l_off` 1-based
start from `@C:123` syntax. `mseq_record`: `lseek` file offset to re-read the
entry; `cont` chunk index (1 = first, `comp_lib9.c:2546`); `libstr` identifier
(kept across chunks, `comp_lib9.c:2550`); `n1tot_p` -> total length of the
original sequence (`comp_lib9.c:2596-2611`); `bline` description.

## `struct getlib_str` (`src/param.h:221-246`) -- library reading state

`lcont`/`ocont` continuation flags, `aa1save` holds the last 150 residues to
prepend to the next chunk (`comp_lib9.c:2505-2508,2644-2646`), `loffset`
running offset, `tot_memK`/`max_memK`/`use_memory` (whole library kept in RAM
for multi-query runs, limit `MAX_MEMK` 16 GB with `BIG_LIB64`+`COMP_THR`,
`defs.h:70-74`), `start_seqr_chain` list of `struct seqr_chain` buffers.

## `struct rstruct` (`src/rstruct.h:5-15`) -- result of one `do_work()`

`score[3]` (ssearch: `[0]` = SW score, `[1]`,`[2]` = 0, `dropgsw2.c:829-830`;
fasta: initn/init1/opt), `valid_stat` (usable for the statistics sample),
`escore`, `comp`/`H` (for `-z 6`), `alg_info` bits (fasta: join done / opt done).

## `struct beststr` (`src/best_stats.h:19-41`) -- a saved hit (one per chunk per strand)

`seq`, `mseq` pointers; `rst`; `frame` (0 forward, 1 reverse strand); `zscore`,
`zscore2`; `repeat_thresh` (raw score threshold for extra HSPs); `a_res` chain of
alignments; `have_ares`; `bbp_link` (another beststr for the same `seq`,
`compacc2.c:2592-2602`). Stored in `bestp_arr[]` up to `MAX_BEST` 60,000
(`defs.h:104`); when full the weakest quarter is dropped and `zbestcut` raised
(`compacc2.c:2578-2585`). Statistics sample `stats[]` is separate, up to
`MAX_STATS` 60,000 (`defs.h:106`, `compacc2.c:2496-2545`).

## `struct a_res_str` (`src/aln_structs.h:35-60`) -- one alignment (HSP)

`sw_score`; `rst` (re-scored); `min0,max0` query bounds, `min1,max1` library
bounds (0-based start, exclusive end, relative to the chunk); `res[]` encoded
alignment (run-length ops); `nres`; `next` -> next non-overlapping alignment in
the same chunk, sorted by score (`wm_align.c:142-146`); `index` position in the
chain (`dropgsw2.c:1033-1037`); `aln_code` (`-m 9c`/CIGAR/BTOP text);
`aln` (`struct a_struct`, `aln_structs.h:6-27`: `d_start0/d_stop0`,
`d_start1/d_stop1` display coordinates, `nident`, `lc`, `ngap_q/ngap_l`).

## `struct f_struct` for Smith-Waterman (`src/dropgsw2.h`)

Per-thread, per-strand state created by `init_work()` (`dropgsw2.c:132-620`):
`ss[]` scalar DP row; `waa_s`/`waa_a` query profiles (`-S` and plain);
`pam2p[2]`; SSE2: `workspace` (3*16*(MAXTST+MAXLIB+32)+256 bytes, ~9.1 MB,
`dropgsw2.c:456-457`), `word_score`/`byte_score` striped profiles, `bias`,
`try_8bit`, `done_8bit`, `done_16bit`; `max_res` alignment buffer size
(`max(3*n0/2, MIN_RES)`, `dropgsw2.c:618`); `kar_p` for `-z 6`.

## Thread buffers (`src/thr_buf_structs.h`, `src/thr_bufs2.h`)

`struct buf_head` holds `buf2_data[]` (pointers to `seq_record`s, frame,
`stats_idx`, `repeat_thresh`), `buf2_res[]` (`rstruct` per entry),
`buf2_ares[]` (alignments), `hdr.buf2_type` (`BUF2_DOWORK`, `BUF2_DOSHUF`,
`BUF2_DOALIGN`). `reader_buf[]`/`worker_buf[]` ring indices are the only shared
mutable state (`thr_bufs2.h:12-27`).
{% endraw %}
