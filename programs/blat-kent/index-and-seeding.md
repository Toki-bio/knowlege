---
title: "blat-kent: Index, tiles, repMatch, hit lookup, clumping"
type: reference
generated: true
---

> Part of the [blat-kent](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/blat-kent/references/index-and-seeding.md)

{% raw %}
# Index, tiles, repMatch, hit lookup, clumping

Paths relative to `src/`. Verified at `4af84d3`. Own calculations are marked `notes`.

## Parameters and defaults

| Param | blat default | gfServer default | Allowed | Where |
|---|---|---|---|---|
| `-tileSize` (K) | 11 DNA; 5 protein-like | 11; 4 with `-trans` | DNA 6..18, protein 3..8 | `blat/blat.c:33`, `:706`, `gfServer/gfServer.c:62`, `:1494`; `jkOwnLib/genoFind.c:617-630` |
| `-stepSize` | = tileSize | = tileSize | any int (no check) | `blat.c:715`, `gfServer.c:1500`, `genoFind.c:642-643` |
| `-minMatch` | 2; 1 protein-like | 2; 3 `-trans` | >= 0 checked as `< 0` only (`blat.c:725`) | `blat.c:35`, `:707`, `gfServer.c:61`, `:1495` |
| `-maxGap` | 2; 0 protein-like | 2; 0 `-trans` | <= 100 | `blat.c:37`, `:709`, `:727-728` |
| `-repMatch` | `gfDefaultRepMatch` (1024 for 11/11) | same | any | `blat.c:734-737`, `gfServer.c:1501-1504`, `genoFind.c:2716-2769` |
| `-oneOff` | off; presence enables | not available | blat only | `blat.c:717`, `genoFind.c:1050-1051`, `:1254-1255` |
| `-minScore` | 30 | n/a (gfClient 30) | | `blat.c:36`, `gfClient/gfClient.c:31` |
| `-minIdentity` | 90 DNA; 25 protein-like | n/a (gfClient 90/25) | percent, float | `blat.c:55`, `:693`, `gfClient.c:32`, `:199-200` |
| `-ooc` | none | none | flat index only | `genoFind.c:684-691` |
| `-noSimpRepMask` | off | off | | `genoFind.c:689-690`, `:1126` |

Usage text says tileSize "Usually between 8 and 12" (`blat.c:96-98`) and repMatch "Typically 256 for
tileSize 12, 1024 for tile size 11, 4096 for tile size 10 ... when stepSize is halved, repMatch is doubled"
(`blat.c:113-118`).

### `gfDefaultRepMatch(tileSize, stepSize, protTiles)` (`genoFind.c:2716-2769`)

DNA table: K=18:2, 17:4, 16:8, 15:16, 14:32, 13:128, 12:256, 11:1024, 10:4096, 9:16384, 8:65536,
7:262144, 6:1048576; then `repMatch = repMatch * tileSize / stepSize` (`:2766-2767`). Protein: K=3:600000,
4:30000, 5:1500, 6:75, 7:10. Examples: (11,11)->1024; (11,5)->2252; (8,8)->65536; (8,4)->131072;
(12,12)->256; (5,5,prot)->1500. gfServer `-trans` sets `repMatch = gfPepMaxTileUse 30000` but it is
overwritten by `gfDefaultRepMatch(4,4,TRUE)=30000` anyway (`gfServer.c:1497`, `:1504`).

## Tile encoding

- DNA tile = 2 bits per base, `a,c,g,t` only; any other char (N, IUPAC, upper case when masking turned it
  to N) returns -1 and the window is skipped (`gfDnaTile`, `genoFind.c:713-726`; lookup table `:695-711`).
- Peptide tile = base-20; `aaVal` < 0 or > 19 invalid, so 'X' (index 20, `lib/dnautil.c:1212`) and 'Z'
  (stop codon from `translateSeq(...,stop=FALSE)`, `lib/dnaseq.c:137-143`; not in `aminoAcidTable`, so
  `aaVal['Z']=-1`, `dnautil.c:1233-1239`) invalidate the window (`gfPepTile`, `genoFind.c:728-742`).
- Flat index when K <= 12 (DNA) / 5 (pep): `tileSpaceSize = 4^K` or `20^K` (`genoFind.c:645-663`).
- Segmented index when K > 12 (DNA): head = 12-mer hash, tail = `K-12` bases stored as 16-bit value in a
  packed `endList` entry of (tail, offset) (`:656-659`, `:133-181`, `gfAddLargeSeq :933-961`,
  `gfSegmentedFindHits :2172-2232`). **`maxPat = BIGNUM`: no repMatch filtering, no `-ooc`, no PCR**
  (`:670-675`, `:686-687`, `:2703-2704`).

## Index build (two passes) (`genoFind.c:1360-1394` in-memory; `:1028-1119` nib/2bit)

1. `gfCountSeq` (`:744-769`): for `i = 0; i <= size-K; i += stepSize` compute tile of `seq[i..i+K)`; if
   `listSizes[tile] < maxPat` increment. Counts saturate at `maxPat`. Only the target grid positions
   (multiples of stepSize per sequence) are ever indexed.
2. `gfAllocLists` (`:841-878`): tiles with `listSizes < maxPat` get space (`count += listSizes`), the rest
   are counted as `overusedCount` and get none. One `needHugeMem(count * sizeof(gfOffset))` block (`:868`).
3. `gfZeroNonOverused` (`:1010-1026`) resets kept tiles' counts to 0 so pass 2 can fill them.
4. `gfAddSeq` (`:904-931`): same grid; `lists[tile][listSizes[tile]++] = offset` for kept tiles (`:923-926`).
   Offsets are concatenated-genome coordinates; `gfSeqSource.start/end` map back (`:1383-1387`).
5. `gfZeroOverused` (`:992-1008`) sets dropped tiles' `listSizes` to 0, so lookups see an empty list.
6. `-ooc`: `oocMaskCounts` pre-sets `listSizes[tile] = maxPat` for every tile in the `.ooc` file so step 1
   never increments them and step 2 drops them (`ooc.c:9-53`, `genoFind.c:688`). Only with an `.ooc` file
   are the 16 period-1/2 simple repeats (`aaaa...`, `acac...`) also masked for DNA (`oocMaskSimpleRepeats`,
   `ooc.c:55-76`, `genoFind.c:689-690`). Without `-ooc` there is no simple-repeat masking of a DNA index.
7. `-makeOoc` (`gfMakeOoc`, `:2469-2549`): counts tiles on **both strands** (`:2508-2531`) with
   `stepSize = tileSize` (`:2474`), writes every tile with count >= `repMatch` (`:2538-2545`). File =
   `oocSig`, tileSize, then 32-bit tile values.
8. `-mask`: for in-memory blat the target is read with masked bases upper-case, indexer sees them as non-ACGT
   (skipped); after indexing the target is unmasked so extension can run through repeats (`blat.c:607-616`,
   usage `:120-126`). For gfServer `-mask` turns masked bases into N before indexing (`genoFind.c:1155-1186`).

### Memory footprint (flat DNA index) `notes`, formula from the allocations above

```
fixed  = tileSpace * (sizeof(gfOffset*) + sizeof(bits32))   # lists ptr array + listSizes  (genoFind.c:678,680)
       = 4^K * 12 bytes on a 64-bit host
entries= (sum over sequences of floor((len-K)/step)+1, minus dropped) * sizeof(gfOffset)   # 4 B (8 B for Huge)
```

| K | step | tileSpace | fixed | tiles for 3 Gbp | entry bytes |
|---|---|---|---|---|---|
| 11 | 11 | 4,194,304 | 50 MB | 272.7 M | 1.09 GB |
| 11 | 5 | 4,194,304 | 50 MB | 600 M | 2.4 GB |
| 11 | 1 | 4,194,304 | 50 MB | 3000 M | 12 GB |
| 12 | 12 | 16,777,216 | 201 MB | 250 M | 1.0 GB |
| 10 | 10 | 1,048,576 | 12.6 MB | 300 M | 1.2 GB |
| 8 | 8 | 65,536 | 0.8 MB | 375 M | 1.5 GB |
| 8 | 4 | 65,536 | 0.8 MB | 750 M | 3.0 GB |

Plus the target itself in RAM for standalone blat (`gfClientSeqList`, 1 byte/base), plus `gfHit`s per
query (lm pool). 32-bit `gfOffset` caps the total at 4 Gbp (`genoFind.c:796-815`, `:1075-1076`); use
`blatHuge`/`gfServerHuge` beyond that (`blat/makefile`, `inc/genoFind.h:68-76`). Raising repMatch keeps
more entries, so memory grows toward the "no dropping" figure above.

## Hit lookup (`genoFind.c:2346-2388` dispatch)

- Flat DNA, no oneOff: `gfFastFindDnaHits` (`:1964-2020`): rolling 2-bit hash over **every query position**
  (`:1988-1993`, uses `ntValNoN` so N maps to a base value rather than skipping); for each position with
  `listSizes[bits] != 0` and not query-masked (`:1998`), every target offset in the list becomes a `gfHit`
  with `diagonal = tStart + qSize - qStart` (`:2001-2013`). Cost ~ sum of list lengths touched.
- `-oneOff` (`gfStraightFindNearHits`, `:2074-2170`): for each query position and each of the K positions,
  substitute each of the other 3 bases and look up (`:2113-2166`), i.e. 3K+1 = 34 lookups per position for
  K=11. Hits are not marked as near-matches; they enter clumping like exact hits. In-memory index only.
- Segmented (`:2172-2232`) walks the 12-mer head list and compares the stored tail.
- `gfFindHitsInRegion` (`:2411-2429`) restricts to one target and `[tMin,tMax)`; used by `-fine`.

## Clumping (`clumpHits`, `genoFind.c:1877-1961`; `clumpNear :1816-1875`; `targetClump :1752-1812`)

1. Hits are distributed into 64 kb buckets by `tStart` (`:1885-1904`).
2. Each bucket is sorted by diagonal (`gfHitSortDiagonal :1601-1628`). Walk sorted hits: a clump continues
   while `hit->diagonal - lastHit->diagonal <= maxGap` (`:1924-1925`), so hits may differ by up to maxGap
   bases of net indel and be arbitrarily far apart along the diagonal within the bucket.
3. Clumps that reach within 300 bp (`gfNearEnough`, `:1814`) of the bucket's end are pushed to the next
   bucket to allow merging (`:1932-1937`); otherwise kept if `clumpSize >= minMatch` (`:1938-1947`).
4. `clumpNear`: within each clump, sort by `tStart`, split where consecutive hits are > 300 bp apart on the
   target (`:1838`; 100 for peptides, `:1824`), keep pieces with >= minMatch hits (`:1840`, `:1863`).
5. `targetClump`: assign `gfSeqSource`; split clumps spanning two sequences (`:1757-1811`).
6. `findClumpBounds`: `qEnd/tEnd = max start + tileSize` (`:1729-1750`). Compute `queryCoverage` (union of
   tile intervals on query, `:1678-1716`) and sort clumps by descending coverage (`:1954`).

`minMatch=1` therefore means any single surviving tile hit becomes a clump (protein-like default).
Note `-maxGap` applies only between consecutive hits of a clump ("Only relevant for minMatch > 1", `blat.c:110`).

## From clump to alignment seeds

- Default DNA: the clump is only a region; `alignComponents` re-aligns query[qStart-250..qEnd+250] vs target
  range (+-500 via `gfiExpandRange`) with `ffFind(ffCdna)` (`gfBlatLib.c:414-462`, `:653`). The index hits
  are **not** reused as seeds here; fuzzyFind finds its own exact "tiles" (`fuzzyFind.c:1159-1232`) and bails
  if the exact-tile alignment is improbable (`evalExactAli > 0.1`, `:1222-1229`).
- `-fine` / `-fastMap` / protein: hits are merged along the diagonal and extended (`clumpToHspRange`
  `gfBlatLib.c:729-833` with `extendHitRight/Left` stopping after `maxDown` 10 (2 for fastMap) points of
  score loss using `dnaScore2`/`aaScore2`; `clumpToExactRange` exact-only for `-fine`, `ffSeedExtend.c:205-289`).

## Seeding probability `notes` (own model, not from the code)

Assume substitutions only, independent per base, identity p, alignment length L on the target grid with
step s: number of grid tiles n = floor((L-K)/s)+1, each matches exactly with probability p^K; a clump needs
>= minMatch hits (ignoring the maxGap/diagonal constraint, which indels tighten further).

P(>= 2 exact 11-mer hits), step 11:

| L | 70% | 75% | 80% | 85% | 90% | 95% |
|---|---|---|---|---|---|---|
| 100 | 0.01 | 0.05 | 0.18 | 0.46 | 0.83 | 0.99 |
| 300 | 0.10 | 0.32 | 0.69 | 0.95 | 1.00 | 1.00 |
| 500 | 0.22 | 0.57 | 0.91 | 1.00 | 1.00 | 1.00 |

Same, L=300, other settings: K=11 s=5: .32/.71/.96/1/1/1; K=10 s=10: .21/.51/.85/.99/1/1; K=9 s=9:
.39/.72/.95/1/1/1; K=8 s=8: .64/.90/.99/1/1/1; K=8 s=4: .93/1/1/1/1/1; K=12 s=12: .05/.19/.52/.89/1/1.
With `-oneOff` (p_tile = p^K + K p^(K-1)(1-p)), K=11 s=11 L=300: .83/.98/1/1/1/1 (vs .10/.32/.69/.95/1/1).
With minMatch=1, K=11 s=11: L=300: .42/.69/.91/.99/1/1.

So the default seed is tuned for >= 90% identity (BLAT's stated mRNA/genome target) and loses most 75%
alignments of 300 bp at the seeding step alone; see `scoring-and-filters.md` for the later filters that
remove the rest, and `repeat-search-notes.md` for how to retune.
{% endraw %}
