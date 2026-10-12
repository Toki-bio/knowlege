# Build, run, test

Paths relative to `src/`. Nothing here was executed in the digest session (read-only sparse checkout);
commands are reconstructed from the makefiles.

## Build

```sh
cd src
# jkOwnLib and blat need lib/jkweb.a (lib/makefile) and inc/common.mk (MACHTYPE := uname -m, COPT=-O3 -g)
make -C lib            # builds lib/$(MACHTYPE)/jkweb.a
make -C jkOwnLib       # builds lib/$(MACHTYPE)/jkOwnLib.a and jkOwnLibHuge.a (-DGFSERVER_HUGE)  (jkOwnLib/makefile)
make -C blat           # build_both: blat and blatHuge (blat/makefile: IDXSIZE=huge, A=blatHuge, symlink blatHuge.c)
make -C gfServer; make -C gfClient
# top level: `make blatSuite` (src/makefile:89) builds the BLAT directories; binaries go to ${BINDIR}=~/bin/${MACHTYPE} (inc/common.mk:467)
```
`inc/userApp.mk` is the per-program rule set (`A = blat`, `preMyLibs += ../lib/$(MACHTYPE)/jkOwnLib.a`,
`blat/makefile`). Whether `lib/` builds without MySQL/other optional deps in this sparse tree was not checked.

Compile-time switches: `-DGFSERVER_HUGE` (64-bit `gfOffset`, `inc/genoFind.h:68-76`); `DEBUG_CLUMP` /
`DEBUG_HITS` in `jkOwnLib/genoFind.c:27-28` dump hits and clumps to stdout; `DEBUG` in `blat.c:643-653`
hard-codes a command line; `dumpBunList` under `DEBUG` in `gfBlatLib.c:1552-1554`.

## Minimal run that exercises the main path

```sh
blat target.fa query.fa out.psl                               # dna/dna, in-memory index
blat -t=dnax -q=dnax target.fa query.fa out.psl               # 6-frame translated both sides
blat -makeOoc=11.ooc -repMatch=1024 genome.2bit /dev/null x   # over-occurring tile file (both strands)
blat -ooc=11.ooc genome.2bit query.fa out.psl
gfServer start localhost 17779 genome.2bit                    # ~minutes to index; prints "Server ready for queries!"
gfServer index genome.untrans.gfidx genome.2bit               # precomputed index for -indexFile / dynserver
gfServer start localhost 17779 -indexFile=genome.untrans.gfidx genome.2bit
gfClient localhost 17779 /path/to/dir query.fa out.psl        # dir must contain genome.2bit
gfServer direct query.fa genome.2bit                          # clumps only, no alignment (gfServer.c:255)
```
Repeat-search oriented example (see `repeat-search-notes.md` for rationale):
```sh
blat -tileSize=8 -stepSize=4 -minMatch=2 -oneOff=1 -repMatch=5000000 -minIdentity=60 -minScore=40 \
     -maxIntron=50 -noHead genome_chunks.fa sine_consensus.fa out.psl
```

## Test suite

- `blat/test/makefile`: `make test` runs `tBasic` (hCrea mRNA/genome in `basic/`: dna, `-prot`, `-t=dnax -q=prot`,
  `-t=dnax -q=rnax`, `-fine`, diff vs `basic/ref*.psl`), `tThrowback` (`throwback/`, `v29skips/`, uses `pslCheck`
  from `src/utils`, not in this checkout), `tIntronMax` (`intron50k/` with `-minScore=190` and `-maxIntron=40000/5000`,
  good worked example of `cutAtBigIntrons`), repeated with `IDXSIZE=64` for `blatHuge`. `hugeTests` need a
  /hive 40 Gb assembly.
- `gfServer/tests/makefile`: `direct*Test`, `prot*`, `trans*`, `pcrTest`, `indexFile*Test`, `dynServer*Test`
  using `input/*.nib|.2bit|.fa` and `expected/*.out|*.psl`; shell drivers in `tests/bin/`.
- `jkOwnLib/tests/freen` only.

Best worked examples: `blat/test/basic` (every mode, small), `blat/test/small/{target,query}.fa + expected.psl`,
`gfServer/tests/input/repeat.fa` (name suggests repeat handling; content not inspected).

## Debugging tips

- `-verbose=2` prints each query name (`blat.c:287`) and non-monotonic stitch warnings (`supStitch.c:202`).
- `gfServer -debugLog` logs clump/hit counts per query (`gfServer.c:388`) and index parameters (`:224-241`).
- `gfServer direct` prints `qStart-qEnd target tStart-tEnd, hits N` per clump (`genoFind.c:1511-1526`) which is
  the quickest way to see whether seeding (not alignment/filtering) is what loses a hit.
- First breakpoints: `clumpHits` (`genoFind.c:1877`) for seeding; `saveAlignments` (`gfBlatLib.c:492`) for score
  filtering; `savePslx` identity test (`gfOut.c:134`) for identity filtering; `ssStitch` loop (`supStitch.c:788`)
  for the 16-cap.
