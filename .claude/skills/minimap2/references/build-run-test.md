# Build, run, test

All line numbers at commit `3c28777`.

## Build

```sh
git clone https://github.com/lh3/minimap2 && cd minimap2
make                      # x86-64: builds SSE2 + SSE4.1 ksw2 objects with runtime dispatch (Makefile:17-22)
make sse2only=1           # single SSE2 build, no dispatch (Makefile:19-20)
make arm_neon=1 aarch64=1 # ARM via sse2neon (Makefile:11-13, 23-30)
make extra                # also builds `sdust` and `minimap2-lite` (example.c) (Makefile:49-54)
make asan=1 | tsan=1      # sanitizer builds (Makefile:32-40)
```

Requirements: C compiler, zlib, pthreads (`LIBS` `Makefile:8`). Default
`CFLAGS = -g -Wall -O2 -Wc++-compat`, `-DHAVE_KALLOC` (`Makefile:1-2`);
asserts are enabled. Output: `minimap2` binary and `libminimap2.a`
(`Makefile:51`). Prebuilt binaries: README.md:81 (`curl ... v2.31 ...tar.bz2`).
Python binding: `pip install .` / `python setup.py build_ext --inplace`
(`setup.py`, `pyproject.toml`).

## Minimal run that exercises the main path

```sh
./minimap2 -a test/MT-human.fa test/MT-orang.fa > test.sam          # README.md:11
./minimap2 -x map-ont -d MT-human-ont.mmi test/MT-human.fa          # README.md:13 (index)
./minimap2 -a MT-human-ont.mmi test/MT-orang.fa > test.sam          # README.md:14 (reuse index)
./minimap2 -c --cs test/MT-human.fa test/MT-orang.fa                # PAF with CIGAR and cs tag
./minimap2 -c test/t-inv.fa test/q-inv.fa                           # inversion rescue toy (tp:A:I lines)
```

Test fixtures in `test/`: `MT-human.fa`, `MT-orang.fa` (mitochondria),
`t-inv.fa`/`q-inv.fa` (inversion), `t2.fa`/`q2.fa`, `x3s-ref.fa`/`x3s-qry.fa`
with `x3s-aln.txt` (expected alignment).

## Test suite

There is **no `make test` target** (targets: `all`, `extra`, `minimap2`,
`minimap2-lite`, `sdust`, `clean`, `Makefile:47-104`). Regression testing
upstream is manual: run the toy fixtures above and `misc/paftools.js`
(`paftools.js mapeval`) on simulated reads. For behavioural checks, compare
PAF output before/after a change with the debugging flags below; the `rl`,
`cm`, `s1` tags expose seeding and chaining without alignment.

A quick repeat-search smoke test (no large genome needed):

```sh
# all-vs-all of a small fasta shows seeding/chaining behaviour with -P
./minimap2 -k11 -w1 -f1000000 -n1 -m20 -P -c test/MT-human.fa test/MT-orang.fa | cut -f1-12,13-16
# which query minimizers were thrown away as repetitive?
./minimap2 --dbg-seed-occ -t1 ref.fa sine_consensus.fa 2>&1 >/dev/null | awk '$1=="SF"' | sort -k2,2n
```

## Debugging tips

- `-v3` (default) prints `[M::mm_mapopt_update::...] mid_occ = N`
  (`options.c:81-82`) and index stats (`index.c:132-133`). This single number
  explains most "no hits in repeats" questions.
- `--dbg-seed-occ` (`main.c:264`): `SF qpos occurrence filtered(0/1)` per query
  minimizer found in the index (`seed.c:115-116`). Works with `-t 1`
  interleaving concerns in mind.
- `--print-seeds` (`main.c:220`, forces one thread): `RS rep_len`, every anchor
  (`SD`) and every chain's anchors (`CN`) (`map.c:255-260, 326-330`).
- `--print-chains` (`main.c:255`): chains only.
- `--print-aln-seq` (`main.c:224`): each ksw2 call's parameters, sequences,
  score and CIGAR (`align.c:339-346, 361-367`).
- `--print-qname` (`main.c:219`): per-query thread, timing and arena use
  (`map.c:433-436, 472-473, 366-373`); good for finding the one query that
  takes minutes.
- `--no-kalloc` (`main.c:218`) swaps the arena for libc malloc so
  ASan/valgrind can attribute leaks (`mm_tbuf_init` `map.c:17`,
  `mm_idx_init` `index.c:58`).
- Breakpoints: `mm_collect_matches` (`seed.c:98`) to see what survives `-f`;
  `mg_chain_backtrack` (`lchain.c:27`) to see why chains fail `-n`/`-m`;
  `mm_select_sub` (`hit.c:255`) for `-p`/`-N` pruning; `mm_filter_regs`
  (`hit.c:301`) for post-alignment drops; `mm_align_pair` (`align.c:336`) for
  DP inputs.
- Deterministic output: ties are broken by a hash of query name, length and
  `--seed` (`map.c:246-248`, `hit.c:64-65`); `--no-hash-name` removes the name
  component.
