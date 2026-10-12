# Build, run, test

## Build

```sh
git clone https://github.com/lh3/bwa && cd bwa
make -j            # gcc -g -Wall -O3 -DHAVE_PTHREAD -DUSE_MALLOC_WRAPPERS; links -lm -lz -lpthread (-lrt on Linux)
make asan=1        # AddressSanitizer build (Makefile:24-27)
make bwamem-lite   # example.c against libbwa.a (Makefile:36-37)
```
Only dependency: zlib (`README.md`). SSE2 is used when `__SSE2__` is defined, NEON via `neon_sse.h`, otherwise a scalar fallback (`ksw.c:29-35`). CI (`.github/workflows/ci.yaml`) builds on x64/arm64 Linux and macOS only; there is **no test target and no test suite**.

## Minimal run that exercises the main path

```sh
./bwa index -p ref ref.fa                       # writes ref.{pac,ann,amb,bwt,sa}
./bwa mem -t 4 ref reads.fq > aln.sam            # SE
./bwa mem -t 4 ref r1.fq r2.fq > aln-pe.sam      # PE
./bwa fastmap -l 17 -w 20 ref reads.fq           # SMEMs only
./bwa shm ref && ./bwa mem ref reads.fq          # index from shared memory
```
Repeat/consensus exploratory run (see `repeat-finding-assessment.md` for rationale):
```sh
./bwa fastmap -l 12 -w 1000000 ref sine.fa | less          # which exact anchors exist and their copy numbers
./bwa mem -k 12 -r 1 -s 1000000 -y 1000000 -c 1000000 -D 0 -A1 -B1 -O1 -E1 -L0 -T 25 -a -Y ref sine.fa > sine.sam
```

Manual index pipeline with a custom SA interval:
```sh
./bwa fa2pac ref.fa ref; ./bwa pac2bwtgen ref.pac ref.bwt; ./bwa bwtupdate ref.bwt; ./bwa bwt2sa -i 8 ref.bwt ref.sa
```
(`bwtindex.c:128` warns the `pac2bwt` output is unusable until `bwtupdate`; `fa2pac` packs forward only unless told otherwise — check `bwa_fa2pac` at `bntseq.c:335` for the `-f` flag before relying on this path.)

## Test suite

None. Behavioural references are the man page (`bwa.1`), `NEWS.md` ("BWA produces alignments identical to v0.7.17 and v0.7.18"), and `example.c` for the library API (`bwa_idx_load`, `mem_opt_init`, `mem_align1`, `mem_reg2aln`).

## Debugging tips

- `-v 4` on `bwa mem` prints, per read: SMEMs/chains (`mem_print_chain`, `bwamem.c:260-275`), "fraction of repetitive seeds" (`337`), every extension with band and score (`740-798`), band inference and final global score (`1141,1148`). Output goes to stdout mixed with SAM: redirect SAM elsewhere or use `-o`.
- `-F 1` (`BWA_DBG_QNAME`) enables qname-level debug prints where implemented (`bwa.h:46`, `fastmap.c:195`).
- `-1` disables the I/O pipeline thread (`fastmap.c:161,395`) for deterministic single-threaded tracing; `-K` fixes batch size.
- First breakpoints: `mem_collect_intv` (`bwamem.c:140`) to see seeds, `mem_chain_flt` (`353`) to see what is dropped, `mem_chain2aln` (`658`) for extension, `mem_mark_primary_se_core` (`519`) for secondary assignment, `mem_approx_mapq_se` (`982`) for MAPQ.
- `bwa fastmap` is the fastest way to inspect seeding independently of alignment.
- `example.c` + `make bwamem-lite` gives a 60-line harness to call `mem_align1` on one sequence under a debugger.
