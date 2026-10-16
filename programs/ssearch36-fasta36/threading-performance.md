---
title: "ssearch36-fasta36: Threading, performance, memory"
type: reference
generated: true
---

> Part of the [ssearch36-fasta36](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/ssearch36-fasta36/references/threading-performance.md)

{% raw %}
# Threading, performance, memory

## Threads (`-T`)

- `-T n` -> `fa_max_workers` (`src/doinit.c:432-444`); default is the number of
  CPUs from `sysconf(_SC_NPROCESSORS_CONF)` (`doinit.c:255`). Windows: 2.
  `-DMAX_WORKERS=8` in `make/Makefile.linux64_sse2:32` only seeds the variable
  (`doinit.c:92`) and sizes nothing else.
- Threaded driver is `comp_lib9.c` compiled `-DCOMP_THR -DCOMP_MLIB`
  (`make/Makefile.pcom_t:128-129`); `pthr_subs2.c` implements the buffer
  hand-off (`get_rbuf`/`put_rbuf` for the reader 156-208, `get_wbuf`/`put_wbuf`
  for workers 322-380). lalign36 is single-threaded (`Makefile.pcom_t:50`).

## How work is split

- Unit of work = one `buf_head` buffer of library entries, not a sequence and
  not a query slice. For DNA a buffer is sized for `max_chain_seqs` =
  `DEF_WORKER_BUF/AVE_NT_LEN` = 6,000,000/1200 = 5000 entries, minimum 600
  (`src/comp_lib9.c:821-830`, `defs.h:84-86`), with a 6 MB residue area
  (`max(max_chain_seqs*ave_seq_len, 2*max_tot)`, `comp_lib9.c:2039`). A buffer
  therefore holds ~40 chunks of 150 kb.
- Number of buffers: between 2x and 4x the worker count (`comp_lib9.c:838-841`).
- Each worker runs `buf_do_work()` serially over its buffer, both strands per
  entry (`compacc2.c:2960-3057`). Load balance is coarse: a 3 Gbp genome is
  ~20,000 chunks = ~500 buffers; with 32 threads each gets ~16 buffers.
- The main thread reads, chunks and runs `save_best()`; reading is sequential
  and single-threaded (`next_seqr_chain`, `comp_lib9.c:2500-2662`), using mmap
  for FASTA files when `-DUSE_MMAP` (`Makefile.linux64_sse2:32`, `mmgetaa.c`).
- Alignments for displayed hits are also distributed (`buf_align_seq()`,
  `comp_lib9.c:1554`; `buf_do_align()` `compacc2.c:3060`), as are shuffles.
- The guide reports >40x speed-up on 48 cores with the "preload" driver
  `comp_lib8.c` (`doc/fasta_guide.tex:2118-2121`); the shipped Makefile uses
  `comp_lib9` (`make/Makefile36m.common:9-10`), which also keeps the library in
  memory for multi-query runs.

## Expected speed

- No GCUPS figure exists in source or docs (grep of `doc/`, `src/`, `make/`).
  The man page claims SSE2 "can speed-up Smith-Waterman calculations 10 -
  20-fold" over scalar (`doc/fasta36.1:67-70`); the 8-bit kernel is "roughly
  50% faster" than 16-bit (`dropgsw2.c:819`).
- `notes` (not from the source): Farrar's 2007 paper reports ~3 GCUPS for the
  8-bit SSE2 kernel on a 2 GHz Xeon core; modern cores at ~4-8 GCUPS/core would
  put a 300 bp query x 3 Gbp genome x 2 strands at ~1.8e12 cells, i.e. ~4-8
  core-minutes for the scan. Measure before planning.
- Cost drivers beyond the kernel: (a) every displayed hit triggers a scalar
  `sw_walign()` over the full 150 kb chunk plus recursive sub-alignments
  (`wm_align.c:37,149`); with thousands of SINE-containing chunks this can
  dominate. Reduce with `-E`, `-b`, or smaller `-N` chunks (less scalar DP per
  hit, more entries). (b) `-z 11..16` doubles the scan. (c) `-z 6` runs
  `do_karlin()` per chunk (`dropgsw2.c:832-838`).

## Memory

- Library residue buffers: `max_work_buf` x 6 MB (`comp_lib9.c:2039`), plus the
  whole library retained up to `max_memK` (`MAX_MEMK` 16 GB for
  `BIG_LIB64`+`COMP_THR`, `defs.h:70-74`; `-XM 4G` to change,
  `initfa.c:1398-1407`) so that a second query does not re-read it
  (`comp_lib9.c:1889-1899`). A 3 Gbp genome fits in 16 GB (1 byte/residue).
- Per worker per strand: SSE2 workspace ~9.1 MB (`dropgsw2.c:456`) + profiles
  `3*(n0+32)*(nsq+1)` bytes + scalar `ss[]` row (`dropgsw2.c:172-178`).
- `bestp_arr` 60,000 `beststr` (`defs.h:104`), `stats` 60,000 (`defs.h:106`).
- Alignment buffers `max_res = max(3*n0/2, MIN_RES=2000)` ints per hit
  (`dropgsw2.c:618`, `defs.h:95`).
- `MAXTST` 40,000 caps the query (`defs.h:92`; truncation `comp_lib9.c:710-715`).

## Compile flags that matter (`make/Makefile.linux64_sse2`)

| Flag | Effect | Where |
|---|---|---|
| `-DSW_SSE2` | compile `dropgsw2.c`, `smith_waterman_sse2.c`, `global_sse2.c`, `glocal_sse2.c`, `dropnnw2.c` with SIMD | `make/Makefile.fcom:211,220,232-250`; selected by `DROPGSW_O = $(DROPGSW_SSE_O)` `Makefile.linux64_sse2:58-61` |
| `-DSW_ALTIVEC` | PowerPC path | `Makefile.fcom:214,223,229` |
| (none) | scalar `FLOCAL_ALIGN` build = `ssearch36s` | `Makefile.pcom_t:62-63`, `Makefile.sse_alt:6` |
| `-DCOMP_THR -DCOMP_MLIB` | threaded, library-in-memory driver | `Makefile.pcom_t:129` |
| `-DMAX_WORKERS=8`, `-DTHR_EXIT=pthread_exit`, `-D_REENTRANT` | thread defaults | `Makefile.linux64_sse2:32` |
| `-DBIG_LIB64 -D_LARGEFILE64_SOURCE -DUSE_FSEEKO -DUSE_MMAP` | >2 GB libraries, mmap reading, 16 GB `MAX_MEMK` | `Makefile.linux64_sse2:32`, `defs.h:70-74` |
| `-DSAMP_STATS` | sub-sample statistics beyond 60,000 | `Makefile.linux64_sse2:32`, `compacc2.c:1331` |
| `-DM10_CONS` | `-m 10` consensus lines | `Makefile.linux64_sse2:32` |
| `-O` only | the shipped CFLAGS are `-g -std=c99 -O`; `notes`: add `-O3 -msse2` (x86-64 implies SSE2) yourself | `Makefile.linux64_sse2:16` |
| ARM | `make -f ../make/Makefile.linux64_simde_arm`; same `-DSW_SSE2` sources through SIMDe | `make/Makefile.linux64_simde_arm:29,54-61` |

No AVX2 flag exists; building with `-mavx2` only lets the compiler widen nothing
(the kernel is written with 128-bit `__m128i`).
{% endraw %}
