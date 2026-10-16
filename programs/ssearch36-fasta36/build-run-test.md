---
title: "ssearch36-fasta36: Build, run, test"
type: reference
generated: true
---

> Part of the [ssearch36-fasta36](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/ssearch36-fasta36/references/build-run-test.md)

{% raw %}
# Build, run, test

## Build

```sh
cd fasta36/src
make -f ../make/Makefile.linux64_sse2 all      # x86-64 Linux, SSE2 kernels, threads
# binaries land in ../bin/: fasta36 ssearch36 lalign36 ggsearch36 glsearch36 ...
make -f ../make/Makefile.linux64_simde_arm all # ARM64 (NEON through SIMDe)
make -f ../make/Makefile.os_x86_64 all         # macOS x86; Makefile.os_x_arm64 for Apple silicon
make -f ../make/Makefile.linux64_sse2 ssearch36s  # scalar reference build
```
(`doc/INSTALL:1-16`, `make/README`, `make/Makefile.pcom_t:48-80` for targets.)
Default CFLAGS are `-g -std=c99 -O` plus the define list at
`make/Makefile.linux64_sse2:16,32`; override with `CFLAGS="-O3 ..."` on the make
line (CFLAGS uses `+=`). `CC ?= gcc`. No configure step, no external deps
beyond libm and pthreads (`THR_LIBS = -lpthread`, `Makefile.linux64_sse2:44`).

## Minimal run that exercises the main path

```sh
# SINE consensus vs genome, both strands, all HSPs with E<=0.01, tabular
ssearch36 -n -T 16 -E "0.01 1" -m 8CBl -z 1 sine.fa genome.fa > hits.m8
# same with clean shuffle null (2x time) and forced db size
ssearch36 -n -T 16 -E "0.01 1" -m 8CBl -z 11 -Z 20000 sine.fa genome.fa
# coordinates + encoded alignment, no alignment blocks
ssearch36 -n -m 9c -d 0 -E "0.01 1" sine.fa genome.fa
# fasta36 heuristic, ktup 4 (positional), show E-threshold fractions in header
fasta36 -n -E "0.01 1" -m 8CBl sine.fa genome.fa 4
# all non-intersecting local alignments between two sequences
lalign36 -n -K 200 -E 0.1 -m 11 sine.fa bac.fa
```
`-n` is required for stdin queries and recommended always (composition
guessing, `fasta36.1:106-110,445-447`). Options must precede the file names
(`fasta36.1:81-93`).

## Test suite

`test/test.sh` (bash) runs every program against `$SLIB2/fa_dbs/qfo20.lseg` and
`%RMB`/`%p` FASTLIBS entries, writing `test/results/*`; it needs those local
databases and a `FASTLIBS` file, so it does not run out of the box.
Useful worked examples inside it:
- ssearch36 with `-m 9c -S -z 3` and `-z 2 -Z 100000 -s P250` (`test/test.sh:27-28`);
- pairwise shuffle statistics `ssearch36 -k 1000 -a q.aa lib.aa` (`test.sh:36`);
- DNA fasta36 with ktup 4 positional and `-r "+2/-4"`, `-c O` classic thresholds
  (`test/test2.sh:39-40`);
- `test/test_z.sh` loops `-z 1 2 3 11 12 21 22` for fasta36 and ssearch36.
`seq/mgstm1.nt` and `seq/mgstm1.aa` are small queries you can use against any
local FASTA file. There are no unit tests; `test/results/` holds expected output.

## Debugging tips

- `-D` (`doinit.c:380`) turns on library range checks and keeps temp files;
  compile with `-DDEBUG` for adler32 CRC checks of every chunk
  (`compacc2.c:3000-3020`) and `st_sort` of statistics.
- `-R stats.file` dumps every chunk's raw scores, length, frame, and
  `stats_idx` (`compacc2.c:2483-2490` format), ideal for checking the chunking
  and the statistics sample.
- `-H` prints the score histogram with the fitted expectation
  (`doinit.c:391-395`); `-Xb` shows z-scores.
- Header line `Smith-Waterman (SSE2, Michael Farrar 2006)` confirms the SIMD
  build (`dropgsw2.c:662-670`); `(PGopt)` means scalar.
- `-m 8C` prints `# Fields:` and `# N hits found` (`mshowbest.c:250-298`) so an
  empty result is distinguishable from a crash.
- First breakpoints: `do_work` (`dropgsw2.c:696`) for scoring,
  `next_seqr_chain` (`comp_lib9.c:2500`) for chunking, `save_best`
  (`compacc2.c:2369`) for what is kept, `nsw_malign` (`wm_align.c:149`) for
  HSP recursion, `process_hist` (`scaleswn.c:236`) for statistics.
{% endraw %}
