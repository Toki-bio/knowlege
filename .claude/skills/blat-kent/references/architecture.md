# Architecture

Full commit SHA: `4af84d30f7920cdd666c4e0c1c9935cf9f4a974b` (merge commit dated 2026-10-11 09:58 -0700,
`git log -1`). Remote: https://github.com/ucscGenomeBrowser/kent. The local clone is a
**sparse, blob-less (`filter=blob:none`) checkout** of `src/blat src/jkOwnLib src/inc src/lib
src/gfServer src/gfClient` only. All paths below are relative to `src/`.

## Directory map

| Path | Role | Key files |
|---|---|---|
| `blat/` | Standalone CLI: builds in-memory index of target, aligns all queries | `blat.c` (761 lines), `makefile` (also builds `blatHuge` with `-DGFSERVER_HUGE`), `test/` |
| `gfServer/` | Persistent in-memory index server + `gfServer index` (mmap-able `.gfidx`) + xinetd "dynamic" server | `gfServer.c`, `design.txt`, `README.blat`, `tests/` |
| `gfClient/` | Asks gfServer for clumps, does alignment locally from `.2bit/.nib` | `gfClient.c` |
| `jkOwnLib/` | BLAT core library (`jkOwnLib.a`, `jkOwnLibHuge.a`) | `genoFind.c` (index, hit lookup, clumping), `gfBlatLib.c` (clump -> alignment driver), `gfOut.c` (output formats), `fuzzyFind.c` (local aligner `ffFind`), `ffSeedExtend.c` (`-fine` path), `supStitch.c` (chaining/stitching), `ffAliHelp.c`, `bandExt.c` (banded affine SW), `trans3.c` (6-frame), `ooc.c` (overused tiles), `gfClientLib.c`, `gfInternal.c/.h`, `gfPcrLib.c` (isPcr) |
| `lib/` | General kent library (`jkweb.a`) | `ffScore.c` (scoring), `ffAli.c`, `dnautil.c` (`dnaScoreMatch`, translation), `dnaseq.c` (`translateSeq`), `psl.c`, `blastOut.c` (blast/blast8/blast9), `axt.c`, `chainBlock.c` (kd-tree chainer), `twoBit.c`, `nib.c` |
| `inc/` | Headers | `genoFind.h` (public API + structs), `fuzzyFind.h`, `supStitch.h`, `aliType.h` (`enum gfType`, `enum ffStringency`), `common.mk`, `userApp.mk` |

## Module dependency sketch

```
blat.c ──┬─> gfClientLib.c   (load fasta/nib/2bit, -mask handling)
         ├─> genoFind.c      gfIndexSeq()            index build
         │                   gfFindClumpsWithQmask() hits -> clumps
         ├─> gfBlatLib.c     gfLongDnaInMem()        chop query, clumps -> bundles -> alignments
         │     ├─> fuzzyFind.c   ffFind()            (default) local aligner
         │     ├─> ffSeedExtend.c ffSeedExtInMem()   (-fine) seed/extend + bandExt.c
         │     └─> supStitch.c   ssStitch()          chain blocks (chainBlock.c for >=10 blocks)
         │           └─> lib/ffScore.c, lib/dnautil.c  scoring
         └─> gfOut.c         gfOutputAny()           psl/pslx (own code), axt/maf/blast* via lib/axt.c, lib/blastOut.c

gfServer.c ─> genoFind.c (gfIndexNibsAndTwoBits / gfIndexTransNibsAndTwoBits, genoFindIndexWrite/Load)
gfClient.c ─> gfBlatLib.c gfAlignStrand()/gfAlignTrans()/gfAlignTransTrans() ─> network ─> gfServer
```

## Entry points

| Command / mode | Entry | Notes |
|---|---|---|
| `blat db query out.psl` (dna/dna, prot/prot) | `blat.c:638` `main` -> `blat.c:559` `blat()` -> `blat.c:613` `gfIndexSeq` -> `blat.c:328` `searchOneIndex` -> `blat.c:240` `searchOne` (both strands) -> `blat.c:201` `searchOneStrand` -> `gfBlatLib.c:1485` `gfLongDnaInMem` | in-memory index of whole target; `-fine` sets `band`, `-fastMap` sets `fastMap` (`blat.c:210`) |
| `blat -t=dnax -q=prot` | `blat.c:621-624` -> `blat.c:451` `bigBlat(transQuery=FALSE)` -> `blat.c:430` `tripleSearch` -> `gfBlatLib.c:949` `gfFindAlignAaTrans` | 3 protein indexes per strand, DB rev-comped between passes (`blat.c:492-554`) |
| `blat -t=dnax -q=dnax/rnax` | `blat.c:625-628` -> `bigBlat(transQuery=TRUE)` -> `blat.c:437` `transTripleSearch` -> `gfBlatLib.c:1569` `gfLongTransTransInMem` -> `gfBlatLib.c:1278` `gfTransTransFindBundles` | 3x3 frame combos per strand |
| `blat -makeOoc=N.ooc` | `blat.c:575-581` -> `genoFind.c:2469` `gfMakeOoc` | counts both strands, exits |
| `gfServer start host port files` | `gfServer.c:1544` -> `gfServer.c:689` `startServer` -> `genoFind.c:442` `genoFindIndexBuild` (or `genoFindIndexLoad` with `-indexFile`) | per-connection: `dnaQuery` `gfServer.c:351`, `transQuery` `:391`, `transTransQuery` `:443` |
| `gfServer index out.gfidx genome.2bit` | `gfServer.c:1595` -> `gfServer.c:1154` `buildIndex` -> `genoFind.c:534` `genoFindIndexWrite` | single `.2bit` only (`gfServer.c:1157-1158`) |
| `gfServer direct probe.fa files` | `gfServer.c:255` `genoFindDirect` | prints clumps only, no alignment |
| `gfClient host port seqDir in.fa out.psl` | `gfClient.c:189` -> `gfClient.c:100` -> `gfBlatLib.c:564` `gfAlignStrand` | client re-reads target `.2bit/.nib` from `seqDir` |

## Design decisions the authors made (and why, if stated)

- Index target not query: "sits in about a gigabyte of memory" (`gfServer/README.blat`), non-overlapping tiles
  keep the index at ~G/K entries (`genoFind.c:758`, `918`).
- Overused tiles are dropped, not down-weighted (`genoFind.c:856-866` comment "If pattern is too much used
  it's no good to us, ignore"). Rationale in usage: speed "by a factor of 40 in many cases" with `-ooc` (`blat.c:69-71`).
- `gfHit`/`endList` are packed to save RAM: endList entry is 6 bytes not 8 "or nearly an extra gigabyte of
  RAM" (`inc/genoFind.h:130-141`).
- Hits sorted with a custom in-line merge sort because "this is the bottleneck for the whole procedure"
  (`genoFind.c:1528-1531`).
- Query is chopped into <= 5000 bp pieces (4500 preferred, 250 overlap) and re-stitched
  (`gfBlatLib.c:1492-1494`, `inc/genoFind.h:473`) to bound fuzzyFind memory/time.
- The aligner is cDNA-vs-genome oriented: `ffCdna` stringency tolerates long target gaps (introns) with
  logarithmic penalties (`lib/ffScore.c:40-72`), slides gaps to GT/AG (`ffAliHelp.c:159-183`), looks for
  small exons (`supStitch.c:230-274`, `ffSeedExtend.c:347-420`).
- Dynamic server: index saved as mmap-able file with offsets not pointers, `madvise(MADV_RANDOM|MADV_WILLNEED)`
  (`genoFind.c:560-583`, `gfServer/design.txt`). "All search parameters, such as tileSize, are stored in the
  index file and can not be changed at runtime" (`design.txt`).
