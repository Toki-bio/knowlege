---
title: "blat-kent: Core data structures"
type: reference
generated: true
---

> Part of the [blat-kent](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/blat-kent/references/data-structures.md)

{% raw %}
# Core data structures

Paths relative to `src/`. Verified at `4af84d3`.

## `struct genoFind` (`inc/genoFind.h:125-165`)

The index. One per (strand, frame) for translated searches (`transGf[2][3]`, `genoFind.h:183`).

| Field | Type | Meaning | Invariant / note |
|---|---|---|---|
| `maxPat` | int | repMatch: tiles with count >= maxPat are dropped | `BIGNUM` when `segSize > 0` (`genoFind.c:673`) |
| `minMatch`, `maxGap` | int | clump thresholds | copied from CLI (`genoFind.c:681-682`) |
| `tileSize`, `stepSize` | int | K and grid spacing in target | `stepSize==0` -> `tileSize` (`genoFind.c:642-643`) |
| `tileSpaceSize`, `tileMask` | int | `4^(tileSize-segSize)` for DNA, `20^(tileSize-segSize)` for peptide | `genoFind.c:645-663` |
| `segSize` | int | 0 = flat index; DNA `tileSize-12` if >12, pep `tileSize-5` if >5 | segmented path uses `endLists` |
| `isPep`, `allowOneMismatch`, `noSimpRepMask` | bool | mode flags | `allowOneMismatch` = `-oneOff` |
| `totalSeqSize` | gfOffset | concatenated target length | 32-bit unless `GFSERVER_HUGE` (`genoFind.h:68-76`) |
| `sources`, `sourceCount` | `gfSeqSource*` | target sequences with `[start,end)` in concatenated coords | bsearch in `findSource` (`genoFind.c:1479-1487`) |
| `listSizes` | `bits32*` | per tile: occurrence count (build) then list length (search); 0 after `gfZeroOverused` for dropped tiles | `genoFind.c:992-1008` |
| `lists` | `gfOffset**` | per tile: pointer into `allocated`, ascending target offsets | flat index only |
| `endLists` | `endListPart**` | per 12-mer head: packed (tail, offset) entries | segmented only, `genoFind.c:133-181` |
| `isMapped` | boolean | mmap'ed from `.gfidx` | don't free `listSizes/allocated` (`genoFind.c:118-122`) |

Created by `gfNewEmpty` (`genoFind.c:632`), filled by `gfSmallIndexSeq`/`gfLargeIndexSeq` (in-memory seqs,
`genoFind.c:1360`, `1396`) or `gfIndexNibsAndTwoBits` (files, `genoFind.c:1028`). Freed by `genoFindFree`.

## `struct gfSeqSource` (`inc/genoFind.h:78-86`)

`fileName` (file-based index, "file.2bit:seqName", `genoFind.c:1104-1105`) or `seq` (in-memory), `start/end`
in concatenated coords, optional `maskedBits`.

## `struct gfHit` (`inc/genoFind.h:88-95`)

`qStart`, `tStart` (concatenated genome coords), `diagonal = tStart + qSize - qStart` (`genoFind.c:2010`).
Allocated in a local-memory pool `lm` per query (`genoFind.c:2007`). The `qSize` term makes diagonals
non-negative; two hits with equal diagonal are collinear with no indel between them.

## `struct gfClump` (`inc/genoFind.h:99-114`)

A run of hits on near-identical diagonals: `qStart/qEnd`, `target`, `tStart/tEnd` (concatenated coords for
blat queries; `tEnd = max(hit.tStart)+tileSize`, `genoFind.c:1729-1750`), `hitCount`, `hitList`,
`queryCoverage` (bases of query covered by tiles, `genoFind.c:1678-1716`). Clump list is sorted by
descending `queryCoverage` (`genoFind.c:1954`). Built by `clumpHits` (`genoFind.c:1877`).

## `struct gfRange` (`jkOwnLib/gfInternal.h:4-22`)

"Roughly corresponds to an exon": `qStart/qEnd`, `tName`, `tSeq`, `tStart/tEnd` (target-relative),
`components` (sub-ranges after `gfRangesBundle`, `gfBlatLib.c:370-410`), `hitCount`, `frame`, `t3`.
Converted from clumps by `seqClumpToRangeList` (`gfBlatLib.c:612-636`) or from HSPs by
`clumpToHspRange` (`gfBlatLib.c:729-833`).

## `struct ffAli` (`inc/fuzzyFind.h:49-59`)

One gapless aligned block in a doubly linked list: `nStart/nEnd` (needle = query pointers),
`hStart/hEnd` (haystack = target pointers), `startGood/endGood`. Gaps between blocks are implicit:
`hGap = right->hStart - left->hEnd`, `nGap = right->nStart - left->nEnd` (`lib/ffScore.c:138-139`).
A psl block list is exactly this list (`gfOut.c:162-169`).

## `struct ssBundle` / `ssFfItem` (`inc/supStitch.h:26-53`)

Bundle = all candidate `ffAli` lists for one query against one target sequence (`qSeq`, `genoSeq`,
`isProt`, `t3List`, `avoidFuzzyFindKludge`). `ssStitch` (`supStitch.c:756`) rewrites `ffList` to the
<= `maxToReturn` best chains. Standalone blat merges bundles by `genoSeq->name` before the final stitch
(`gfBlatLib.c:1330-1354`).

## `struct gfOutput` (`inc/genoFind.h:206-252`)

Polymorphic writer: `out()` per alignment, `queryOut()` per query (axt/blast bundles flushed here),
`fileHead()`. `minGood` = minIdentity in parts per thousand (`gfOut.c:444`, set from `minIdentity*10`
at `blat.c:589`). `maskHash` separates `repMatch` vs `match` counts in psl (`gfOut.c:92-99`).

## `struct genoFindIndex` + file headers (`inc/genoFind.h:174-184`, `genoFind.c:196-221`, `462-476`)

Container for one untranslated or six translated `genoFind`s; on-disk `.gfidx` = `genoFindIndexFileHdr`
(magic "genoFind", version "1.0", `indexAddressSize` 32/64) followed by per-index `genoFindFileHdr`
with offsets to sources, listSizes, lists/endLists; mmap'ed whole (`genoFind.c:560-583`).

## Enums (`inc/aliType.h`)

`gfType`: dna=0, rna=1, prot=2, dnax=3, rnax=4 (`:6-14`). `ffStringency`: ffExact=0, ffCdna=1 (long target
gaps ok), ffTight=2, ffLoose=3 (`:22-30`). Constants: `gfMinMatch 2, gfMaxGap 2, gfTileSize 11,
gfMaxTileUse 1024, gfPepMaxTileUse 30000` (`inc/genoFind.h:60-66`); `MAXSINGLEPIECESIZE 5000` (`:473`);
`ffIntronMaxDefault 750000` (`inc/fuzzyFind.h:63`).
{% endraw %}
