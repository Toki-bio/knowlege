---
title: "ssearch36-fasta36: Architecture"
type: reference
generated: true
---

> Part of the [ssearch36-fasta36](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/ssearch36-fasta36/references/architecture.md)

{% raw %}
# Architecture

Full commit SHA: `6d6718aaf949483197dd48ea9baecfb7d0f73a41` (2026-09-15, "fix typos,
better error message, remove URI::encode (redundant)").

## Directory map

| Path | Role | Key files |
|---|---|---|
| `src/` | all C sources, one flat directory (~74k lines) | see below |
| `make/` | per-platform Makefiles; program/link rules shared | `Makefile.linux64_sse2`, `Makefile36m.common`, `Makefile.pcom_t`, `Makefile.fcom`, `Makefile.sse_alt` |
| `doc/` | man page, LaTeX guide, change log | `fasta36.1`, `fasta_guide.tex`, `readme.v36`, `fasta.options` |
| `test/` | shell scripts that run every program against local DBs | `test.sh`, `test_z.sh`, `test2.sh` |
| `seq/` | tiny query sequences used by tests | `mgstm1.aa`, `mgstm1.nt` |
| `scripts/`, `psisearch2/`, `misc/` | annotation/expansion scripts, PSI-Search, lav2ps | not needed for search internals |
| `src/simde/` | vendored SIMD-everywhere headers (x86 intrinsics on ARM) | `simde/x86/sse2.h` |

## Which sources make which program

Every program = shared driver + one `drop*.c` comparison module + an `initfa.c`
compiled with a program define (`make/Makefile.fcom:8-30`). Link lines are in
`make/Makefile.pcom_t:50-80`.

| Program | Comparison module (compile defines) | Display/consensus | Notes |
|---|---|---|---|
| ssearch36 | `dropgsw2.c` + `smith_waterman_sse2.c` (`-DSW_SSE2`, `Makefile.fcom:211,232`) + `wm_align.c` | `cal_cons.c` (calcons_sw) | `DROPGSW_SSE_O`, `Makefile.sse_alt:7` |
| ssearch36s | `dropgsw2.c` scalar (`DROPGSW_NA_O`, `Makefile.sse_alt:6`) | same | `FLOCAL_ALIGN` C loop `dropgsw2.c:875-990` |
| lalign36 | `dropgsw2.c` with `-DLALIGN` + `lsim4.c` (`Makefile.fcom:217-220`) | `cal_cons2.c -DLALIGN -DLCAL_CONS` | non-threaded driver `comp_mlib9.o` (`Makefile.pcom_t:50`) |
| fasta36 | `dropnfa.c` (`DROPNFA_O`, `Makefile.pcom_t:6`) + `wm_align.c` | `cal_cons2.c -DFASTA` | banded alignments by default for DNA |
| ggsearch36 / glsearch36 | `dropnnw2.c` (`-DGLOBAL_GLOBAL` for gg) + `global_sse2.c` / `glocal_sse2.c` | `cal_cons.c` | `-z 0` default (`initfa.c:936-938`) |
| fastx/y36, tfastx/y36 | `dropfx2.c`, `dropfz3.c` | `cal_consf.c` | translated; out of scope |

Driver objects: `comp_mthr9.o` = `comp_lib9.c -DCOMP_THR -DCOMP_MLIB`
(`Makefile.pcom_t:128-129`), `work_thr2.o` (worker thread body),
`pthr_subs2.o` (mutex/condvar buffer exchange), `compacc2_t.o` (= `compacc2.c`,
save_best/buf_do_work/reset_maxn), `showbest.o` (= `mshowbest.c`),
`mshowalign2_t.o`, `scale_se.o` (= `scaleswn.c`, statistics), `karlin.o`,
`lgetlib.o`/`lgetaa_m.o` (= `nmgetlib.c`, `mmgetaa.c` readers), `doinit.o`
(generic options), `init_sw_sse.o` (= `initfa.c -DSSEARCH -DSW_SSE2`).

## Module dependency sketch (ssearch36)

```
main()  comp_lib9.c
  |- initenv()        doinit.c  -> f_initenv/f_getopt/parse_ext_opts  initfa.c  (options -> pstruct, mngmsg)
  |- resetp()/last_init() initfa.c  (DNA vs protein defaults, matrix, Karlin lambda/K)
  |- init_getlib_info()/next_seqr_chain()  comp_lib9.c  -> GETLIB = agetlib()  nmgetlib.c  (reads, chunks)
  |- init_thr()  pthr_subs2.c  -> work_thread()  work_thr2.c
  |      |- init_aa0() -> init_work()  dropgsw2.c   (build striped profile per thread, per strand)
  |      |- get_wbuf() -> buf_do_work()  compacc2.c -> do_work()  dropgsw2.c -> smith_waterman_sse2_byte/word()
  |      |- buf_shuf_work() (-z 11..)      buf_do_align() -> do_walign() dropgsw2.c -> nsw_malign() wm_align.c
  |- save_best()  compacc2.c   (one beststr per chunk/strand, stats sample)
  |- process_hist()/scale_scores()/zs_to_E()  scaleswn.c
  |- buf_shuf_seq()  comp_lib9.c   (shuffle fallback / -z 2x)
  |- showbest()  mshowbest.c   (-m 8/9 rows, one per HSP)
  |- showalign() mshowalign2.c -> calc_cons/ calc_code  cal_cons.c, c_dispn.c
```

## Entry points

| Command / mode | Entry | Notes |
|---|---|---|
| any program | `main()` `src/comp_lib9.c:~440` (variables at 452-560; option parsing `initenv()` at `doinit.c:~330`) | one binary per algorithm, chosen at link time |
| `-X?` extended options | `parse_ext_opts()` `src/initfa.c:1346-1440` | `-XM`, `-Xo`, `-Xy`, `-XA`, `-Xg`, `-XB` |
| `-m` formats | `pre_parse_markx()`/`parse_markx()` `src/doinit.c:~760,875-955` | `-m "F8C out.file"` writes extra files |

## Design decisions the authors made (and why, if stated)

- One driver, many `drop*.c` modules with a fixed function interface
  (`init_work`, `do_work`, `do_opt`, `do_walign`, `close_work`, declared in
  `drop_func.h`) so statistics, threading and output are shared.
- Scan and align are separate passes: the scan only needs the score, so the SIMD
  kernel never does traceback; alignments are made only for displayed hits
  (`comp_lib9.c:1524-1556`, `mshowbest.c:414-465`).
- FASTA36 added multiple HSPs per library sequence (previous versions showed only
  the best) by recursive sub-alignment around the best alignment
  (`doc/readme.v36:2-7`, `doc/README.versions:7-19`).
- E-values from the empirical score distribution of the actual library
  (regression on length) rather than from precomputed Karlin-Altschul tables
  (`scaleswn.c:30-40`), with shuffle fallbacks when the library is too small
  (`comp_lib9.c:1077-1096`).
- Long sequences are chunked with overlap so buffer sizes stay bounded
  (`defs.h:91-93`, `doc/fasta_guide.tex:2126-2136`).
{% endraw %}
