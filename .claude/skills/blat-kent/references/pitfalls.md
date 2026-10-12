# Pitfalls and gotchas

Paths relative to `src/`. Verified by reading at `4af84d3` unless marked otherwise; nothing was executed.

- **At most 16 alignments per (query, strand, target sequence) in standalone blat.** `ssAliCount = 16`
  (`jkOwnLib/gfBlatLib.c:22`); `addToBigBundleList` merges all bundles that share `genoSeq->name`
  (`:1330-1354`, called `:1549`); the final `ssStitch(bun, ffCdna, minScore, ssAliCount)` (`:1557`) stops after
  `maxToReturn` chains (`jkOwnLib/supStitch.c:837-841`). A SINE with 100k copies on chr1 yields <= 16 psl lines
  per strand for chr1. Workaround: split the target into many small named sequences (the cap is per name), or
  change `ssAliCount`. gfClient path caps 16 per `gfRangesBundle(ffIntronMax)` cluster instead (`:579`, `:589`)
  but is bounded by the server's `-maxDnaHits` (100 clumps, `gfServer/gfServer.c:68`, `:383-384`).
- **repMatch drops every tile with >= repMatch occurrences on the target grid**, saturating count
  (`jkOwnLib/genoFind.c:762-765`, `:857`, `:1002`). Default 1024 for 11/11 (`:2751-2752`). High-copy repeat
  families lose nearly all seeds. Raise `-repMatch` (memory grows linearly with kept occurrences) or build
  the index only on a subset.
- **`-oneOff=N` is a presence flag**: `oneOff = optionExists("oneOff")` (`blat/blat.c:717`) despite the option
  table declaring `OPTION_INT` (`:172`) and usage saying "If set to 1" (`:100-101`). `-oneOff=0` turns it on.
  Not supported by gfServer at all (`genoFind.c:1050-1051`, `:1254-1255`; the server never sets it, `gfServer.c:65`).
- **Simple-repeat tile masking for DNA happens only when an `-ooc` file is given** (`genoFind.c:684-691`).
  Without `-ooc`, `aaaaaaaaaaa` etc. are indexed like any tile (unless dropped by repMatch). Protein indexes
  always mask period-1/2 repeats (`:1121-1153`, `:1371-1372`).
- **tileSize > 12 disables repMatch filtering and `-ooc`** (segmented index, `maxPat = BIGNUM`, `:670-675`,
  `:686-687`). Memory per entry rises to 6 bytes (`:140`).
- **Default `-maxIntron=750000` lets the stitcher chain distant repeat fragments into one psl line**; cutting
  happens only when `-maxIntron` is non-default (`supStitch.c:815-820`), and the ffCdna gap penalty for a
  500 kb target gap is only ~43 (`lib/ffScore.c:40-72`). For repeat searches pass a small `-maxIntron`.
- **minIdentity is on block bases with target gaps free (ffCdna)**, numerator subtracts 2 per query gap count
  (`jkOwnLib/gfOut.c:131-135`). It is applied in the psl writer, after alignment; `-minScore` is applied before
  (`gfBlatLib.c:508`) and also capped at 20 inside the stitcher (`supStitch.c:775-776`), so `-minScore` above
  20 does not prune chaining, only output.
- **Default aligner is fuzzyFind, not seed-extend.** Index hits only define the region; `ffFind` re-seeds with
  its own improbable-oligo tiles and aborts when `evalExactAli > 0.1` (`jkOwnLib/fuzzyFind.c:1222-1229`).
  `-fine` switches to real seed extension of the index hits (`ffSeedExtend.c:1228-1261`).
- **Query is split into <= 5000 bp pieces** (4500 + 250 overlap, `gfBlatLib.c:1492-1494`); `-fastMap` refuses
  queries > 5000 (`blat.c:205-209`). Queries > 5 Mb only warn (`blat.c:29`, `:300-304`).
- **N in query**: `gfFastFindDnaHits` uses `ntValNoN` (`genoFind.c:1984`, `:1990`) so N does not break the
  rolling hash (it maps to some base value); in the target, N (`gfDnaTile`) invalidates the tile (`:721-722`).
- **Stop codons in `-t=dnax`/`-q=dnax`** become 'Z' (`lib/dnaseq.c:137-143`) which is not an indexable amino
  acid (`genoFind.c:737-738`, `lib/dnautil.c:1233-1239`): every 5-aa window overlapping a stop is unseeded.
  Non-coding sequence has a stop roughly every 21 codons per frame.
- **gfServer and blat compute `repMatch` the same way** since 34x5 (`blat/version.doc`), but `gfServer -trans`
  default tileSize is 4 and minMatch 3 (`gfServer.c:1491-1498`) vs blat's 5 and 1 (`blat.c:704-710`): results
  differ between `blat -t=dnax` and `gfServer -trans` + `gfClient`.
- **Target names**: file-based index names targets `file.2bit:seq` (`genoFind.c:1104-1105`); `gfClient` strips
  to the sequence name unless `includeTargetFile` (`gfBlatLib.c:550-562`).
- **Index requires a single contiguous allocation** for all entries (`needHugeMem(count*4)`, `genoFind.c:868`);
  the `int` return of `gfAllocLists` (`:841`) can overflow above 2^31 entries though the value is unused (`:1078`).
- **blast8/9 e-values use a hard-coded 3.0e9 database size** (`lib/blastOut.c:115-121`), and gfClient hard-codes
  `databaseLetters = 3.0e9` (`gfClient/gfClient.c:116-117`).

## Defaults that surprise

| Parameter | Default | Why it matters |
|---|---|---|
| `-minIdentity` | 90 (DNA), 25 (prot/dnax) (`blat.c:55`, `:693`) | removes all 70-85 % hits after they were computed |
| `-repMatch` | 1024 (11/11), scaled `tileSize/stepSize` (`genoFind.c:2766-2767`) | kills seeds of high-copy repeats |
| `-maxIntron` | 750000 (`inc/fuzzyFind.h:63`) | chains far-apart fragments; `cutAtBigIntrons` off at default |
| `-stepSize` | = tileSize (`blat.c:715`) | target grid only; query is scanned every base |
| `-minMatch` | 2 DNA / 1 protein-like (`blat.c:35`, `:707`) | single-hit clumps for protein modes |
| `-maxGap` | 2 (`blat.c:37`) | net indel tolerance between consecutive seeds in a clump |
| `ssAliCount` | 16, not a flag (`gfBlatLib.c:22`) | per-target-sequence output cap |
| gfServer `-maxDnaHits` | 100 (`gfServer.c:68`) | clumps per query over the wire |
| gfServer `-maxNtSize` | 40000 (`gfServer.c:58`) | longer queries rejected server-side |
| `-trimA` | on for `-q=rna/rnax` (`blat.c:684-685`) | polyA removed from RNA queries |

## Known bugs / open upstream issues (as of pinned commit)

- Suspected (not run): `ENDLIST_ENTRY_NUM_PARTS` is 3 under `GFSERVER_HUGE` and 5 otherwise
  (`genoFind.c:134-138`), but `endListSetEntryValues` writes `entry[0..3]` (4 parts) for HUGE and
  `entry[0..2]` (3 parts) otherwise (`:151-165`); `inc/genoFind.h:138-141` documents 48-bit (3x16) for 32-bit and
  80-bit (5x16) for 64-bit. The `#if` branches look swapped: `blatHuge`/`gfServerHuge` with DNA tileSize > 12
  (or protein > 5) would write past each 6-byte entry. Non-huge builds merely waste 4 bytes/entry.
- `clumpToHspRange` duplicate condition `qs != lastQs || ... || qs != lastQs` (`gfBlatLib.c:791`) never
  compares `te`; harmless.
- Historical fixes worth knowing from `blat/version.doc`: segfault with `maxIntron <= 11` and `tileSize <= 11`
  (34x11); `-fine -minIdentity=0 -minScore=0` segfault (34x12); `-minScore` filter broken when over half query
  size (34x3).
