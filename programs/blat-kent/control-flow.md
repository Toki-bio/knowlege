---
title: "blat-kent: Control flow"
type: reference
generated: true
---

> Part of the [blat-kent](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/blat-kent/references/control-flow.md)

{% raw %}
# Control flow

Paths relative to `src/`. Verified at `4af84d3`.

## `blat database query out.psl` (DNA query vs DNA target, default)

1. `main` (`blat/blat.c:638`): `optionInit`; `-t/-q` -> `gfType` (`:661-699`); protein-like defaults
   (`:704-710`); `tileSize`, `stepSize=optionInt("stepSize", tileSize)` (`:714-715`), `minMatch`,
   `oneOff=optionExists` (`:717`), `fastMap`, `minScore`, `maxGap` (<=100, `:727`), `minIdentity`;
   `gfCheckTileSize` (`:724`); `repMatch` = `-repMatch` or `gfDefaultRepMatch(tileSize, stepSize, protLike)`
   (`:734-737`); `setFfIntronMax(-maxIntron)` (`:754`); `setFfExtendThroughN` (`:755`).
2. `blat()` (`:559`): `gfClientFileArray` expands db/query args (single fa/nib/2bit or list file,
   `jkOwnLib/gfClientLib.c:16-52`). `-makeOoc` short-circuits to `gfMakeOoc` and exits (`:575-581`).
3. `gfClientSeqList` loads the whole target into memory as `dnaSeq` list, lower case = unmasked
   (`:583`, `gfClientLib.c:160-208`); counts `databaseSeqCount/Letters` for blast output (`:585-587`).
4. `gfOutputAny(outputFormat, minIdentity*10, ...)` (`:589-590`, `jkOwnLib/gfOut.c:535-577`).
5. If `-mask` unset, target is fully unmasked before indexing (`:611-612`); `gfIndexSeq(dbSeqList, minMatch,
   maxGap, tileSize, repMatch, ooc, isProt, oneOff, FALSE, stepSize, noSimpRepMask)` (`:613-614`) ->
   `gfNewEmpty` + `gfSmallIndexSeq` (`jkOwnLib/genoFind.c:1444-1465`, `1360-1394`): count tiles
   (`gfCountSeq :744`), allocate lists skipping overused (`gfAllocLists :841`), zero non-overused counts
   (`:1010`), add offsets (`gfAddSeq :904`), zero overused counts (`:992`). Target then unmasked again
   (`:615-616`) because extension always runs on unmasked sequence (comment `:607-610`).
6. `searchOneIndex` (`:328`): writes header (`gfOutputHead :337`), iterates query records from nib/2bit/fa
   (`:341-396`); each via `searchOneMaskTrim` (`:308`): `maskQuerySeq` (case handling, `-qMask`, 5 Mb warning,
   `:280-306`), `trimSeq` (polyA/T, `:259-277`), u->t for rna (`:320-321`).
7. `searchOne` (`:240`): `searchOneStrand(seq, isRc=FALSE)`, `reverseComplement`, `searchOneStrand(isRc=TRUE)`,
   `reverseComplement` back, then `gfOutputQuery` (`:251-256`).
8. `searchOneStrand` (`:201`): `-fastMap` size check (<= 5000, `:205-209`); `gfLongDnaInMem(seq, gf, isRc,
   minScore, qMaskBits, gvo, fastMap, band=optionExists("fine"))` (`:210`).
9. `gfLongDnaInMem` (`jkOwnLib/gfBlatLib.c:1485-1566`): chop query into pieces (whole if <= 5000, else 4500
   with 250 overlap, `:1492-1525`). Per piece:
   - default: `gfFindClumpsWithQmask` (`:1537`) -> `gfClumpsToBundles` (`:1544`);
   - `-fastMap`: `fastMapClumpsToBundles` (`:1540`);
   - `-fine`: `ffSeedExtInMem` (`:1533`).
   `addToBigBundleList` merges bundles by target name (`:1549`, `:1330-1354`).
10. `gfFindClumpsWithQmask` (`genoFind.c:2390-2409`): `gfFindHitsWithQmask` dispatch (`:2346-2388`) ->
    `gfFastFindDnaHits` (flat DNA index, no oneOff, `:1964-2020`) -> `clumpHits` (`:1877-1961`) ->
    `clumpNear` (`:1816-1875`) -> `targetClump` (`:1752-1812`) -> coverage sort (`:1953-1954`).
11. `gfClumpsToBundles` (`gfBlatLib.c:638-667`): clumps -> ranges (`seqClumpToRangeList :612`), sort by
    target, `gfRangesBundle(rangeList, 2000)` (`:649`) groups ranges < 2000 bp apart on same target;
    `gfiExpandRange(+-500)` (`:653`, `gfInternal.c:23-54`); per bundle `alignComponents` runs
    `ffFind(..., ffCdna)` on each component +-250 bp (+500 at ends) (`:414-462`); `ssStitch(bun, ffCdna,
    minScore, 16)` (`:661`).
12. Final per-target-sequence: `ssStitch(bun, ffCdna, minScore, 16)` (`:1557`), `refineSmallExonsInBundle`
    (`:1558-1559`, not for fastMap/fine), `saveAlignments` (`:1560-1561`, `:492-514`): `scoreAli` per
    alignment, emit if `score >= minScore` via `out->out` = `pslOut` (`gfOut.c:196-209`) -> `savePslx`
    (`:36-194`) which applies the identity filter (`:131-135`).

## `ssStitch` (chaining) (`jkOwnLib/supStitch.c:756-845`)

1. `minScore` capped at 20 (`:775-776`). Concatenate all `ffAli` lists, sort needle-first, `ffMergeClose`
   (merge same-diagonal blocks that overlap or abut within 3 bp, `ffAliHelp.c:115-156`).
2. Loop: `ssFindBest` (`:703-719`): < 10 blocks -> DAG dynamic program (`ssGraphMake :328`,
   `ssDynamicProgram :428`), else kd-tree `chainBlocks` (`:602-701`, `lib/chainBlock.c`) with
   `ssConnectCost` = overlap adjustment + `ffCalcGapPenalty` (`:563-599`, `:543-552`).
3. Best chain: `ffMergeNeedleAlis`, `ffMergeHayOverlaps`, `forceMonotonic` (`:210-228`), first time only
   `smallMiddleExons` (re-runs `ffFind` in gaps >= 3 bp both sides, `:230-274`, skipped when
   `avoidFuzzyFindKludge`), `cutAtBigIntrons` only if `-maxIntron` != default (`:815-820`, `:721-754`),
   `ffSlideIntrons`.
4. Keep if `score >= minScore`, else stop; stop after `maxToReturn` chains (`:824-841`). Leftover blocks
   feed the next iteration.

## `-fine` (`jkOwnLib/ffSeedExtend.c:1228-1261`)

1. Same clumps (`:1238`). `clumpToExactRange`: merge diagonally adjacent hits, extend exactly both ways
   (`:205-289`). `gfRangesBundle(rangeList, ffIntronMax)` (`:1242`).
2. Per bundle: `ssStitch(bun, ffCdna, 16, 10)` (`:1254`), then `refineBundle` (`:1199-1225`):
   `expandGapless` (+1/-3 per base, stop after drop of 9, `:500-530`, `:52-113`) and restitch;
   `scanIndexForSmallExons` (re-query index for query gaps 11..200 bp, `:347-420`); `bandedExtend`
   (banded affine SW, `axtScoreSchemeRnaDefault` = match 100, mismatch -200, gapOpen/extend 300,
   maxInsert 3, band <= 2*minGap, `:424-496`, `lib/axt.c:493-501`); `scanForSmallerExons`,
   `refineSpliceSites`, `scanForTinyInternal`, `smoothSmallGaps`, `trimFlakyEnds`.

## `-fastMap` (`gfBlatLib.c:860-885`)

Clumps -> HSPs directly (`clumpToHspRange` with `maxDown 2`, `minSpan 20`, `:747-751`), no fuzzyFind;
`gfRangesBundle(rangeList, 256)`; chained by the common final `ssStitch`. Query must be <= 5000 bp (`blat.c:205-209`).

## `-t=dnax` translated target

1. `bigBlat` (`blat.c:451-556`): for `isRc` in {0,1}: `seqListToTrans3List` builds 3 frames per target seq
   (`:403-428`, `trans3.c:9-33`, stops -> 'Z' via `translateSeq(seq,pos,FALSE)`, `lib/dnaseq.c:137-143`);
   `gfIndexSeq` per frame with `isPep=TRUE` (`:504-505`), peptide tiles 20^K, `maskSimplePepRepeat`
   always applied unless `-noSimpRepMask` (`genoFind.c:1121-1153`, `:1371-1372`).
2. `-q=prot`: `gfFindAlignAaTrans` (`gfBlatLib.c:949-995`): clumps per frame -> `clumpToHspRange`
   (aa scoring, `maxDown 10`) -> bundle within `ffIntronMax/3` -> `ssStitch(ffCdna)` -> `saveAlignments`.
3. `-q=dnax/rnax`: `gfLongTransTransInMem` (`:1569-1623`, pieces of 1200 with 270 overlap) ->
   `gfTransTransFindBundles` (`:1278-1328`): query translated 3 ways, 3x3 frame clump sets, HSPs
   untranslated to nucleotide coords (`untranslateRangeList :1154`), bundle within 2000, stitch with
   `ffLoose` for dnax / `ffCdna` for rnax (`:1292`). Then target rev-comped and repeated (`blat.c:550-553`).

## gfServer / gfClient

1. `gfServer start` (`gfServer/gfServer.c:689`): `genoFindIndexBuild` (`:721-722`, oocFile NULL) or
   `genoFindIndexLoad(-indexFile)` (`:727`); accept loop (`:740`), per connection reads one command
   (`:770`): `query|protQuery|transQuery` (`:842-917`, size limit `maxNtSize 40000 / maxAaSize 8000`,
   `:58-59`, `:881`), `pcr`, `status`, `files`, `quit`.
2. `dnaQuery` (`:351-389`): `gfFindClumps`; sends one line per clump `qStart qEnd fileName tStart tEnd hitCount`
   (`:370-372`); stops after `maxDnaHits` (100) clumps (`:361-362`, `:383-384`) or `maxDnaHits/2` per
   sequence with `-perSeqMax` (`:376-382`).
3. `gfClient` (`gfClient/gfClient.c:100-187`): per query, `gfAlignStrand` both strands (`:149-154`) ->
   `gfBlatLib.c:564-596`: `gfQuerySeq` reads clump lines into `gfRange`s (`:136-171`), `gfRangesBundle(ffIntronMax)`
   (`:579`), loads target slice +-500 from `.2bit/.nib` (`gfiExpandAndLoadCached`, `:583`), `alignComponents`
   (`ffFind`), `ssStitch(16)`, `saveAlignments`. Note gfClient hard-codes `databaseLetters 3.0e9` for blast
   e-values (`gfClient.c:116-117`).

## Threading / parallelism model

Single-threaded everywhere. gfServer handles one connection at a time (`gfServer.c:740-...`); the dynamic
server gets parallelism from xinetd spawning processes sharing an mmap'ed index (`gfServer/design.txt`).

## Error handling conventions

`errAbort` aborts the process. fuzzyFind uses `setjmp/longjmp` to recover from allocation failure
(`fuzzyFind.c:1440-1452`). gfServer wraps queries in `errorSafeQuery` (`gfServer.c:575`) so a bad query does
not kill the server; gfClient catches server "Error:" lines as warnings (`gfBlatLib.c:154-158`).
{% endraw %}
