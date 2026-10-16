---
title: "blat-kent — digest"
type: reference
topic: Program digests
generated: true
---

> Type: reference · Claude Code skill `blat-kent` · [source on GitHub](https://github.com/Toki-bio/knowlege/tree/main/.claude/skills/blat-kent)

**When to use:** Internals of BLAT (UCSC kent tree: blat, gfServer, gfClient, jkOwnLib genoFind/fuzzyFind/supStitch). Use when reading, modifying, debugging or explaining BLAT source, its index (tileSize/stepSize/repMatch), seeding and clumping, scoring and -minIdentity/-minScore filters, psl/pslx/blast8 output, gfServer tuning, or when deciding whether and how to use BLAT as a fast pre-filter for diverged repeats (e.g. SINE consensus vs whole genome).

---

{% raw %}
# BLAT (kent) internals

> Type: reference (facts verified against the pinned commit) with marked `notes`.
> Behavioral claims carry `file:line` valid at the pinned SHA. Paths are relative to `src/`.

## Pinned source

| Field | Value |
|---|---|
| Repo | https://github.com/ucscGenomeBrowser/kent |
| Version / tag | BLAT `gfVersion "39x1"` (`inc/genoFind.h:475`); no git tag |
| Commit | `4af84d3` (full SHA in references/architecture.md) |
| Digested on | 2026-10-11 |
| Checkout | **sparse, blob-less shallow clone** at `/home/user/src/kent`: only `src/{blat,jkOwnLib,inc,lib,gfServer,gfClient}`. `src/hg`, `src/utils` (pslReps, pslCDnaFilter, faToTwoBit) are absent. |
| Language / build | C99, GNU make (`inc/common.mk`, `inc/userApp.mk`); `-DGFSERVER_HUGE` builds 64-bit-offset `blatHuge`/`gfServerHuge` |
| License | `blat/LICENSE`, `jkOwnLib/LICENSE`: free for personal/academic/non-profit; commercial use needs agreement with Jim Kent. `lib/` is UC Regents (see `lib/ffScore.c:3-4`). |

## What it does in three sentences

BLAT indexes every non-overlapping K-mer ("tile", default K=11) of the target genome into a
hash from tile value to sorted list of genome offsets, dropping tiles that occur `repMatch`
or more times. A query is scanned at every base; tile hits are bucketed by genome position,
sorted by diagonal, and runs of >=`minMatch` hits within `maxGap` of one diagonal become
"clumps" that are handed to a cDNA-oriented local aligner (fuzzyFind, or seed-extend with
`-fine`) and a chainer (supStitch). Output is filtered by score (`-minScore`) and percent
identity (`-minIdentity`), then written as psl/pslx/axt/maf/blast/blast8/blast9.

## Where to look (read on demand)

| Question | File |
|---|---|
| Code layout, who calls whom, entry points | `references/architecture.md` |
| `genoFind`, `gfHit`, `gfClump`, `gfRange`, `ffAli`, `ssBundle` | `references/data-structures.md` |
| Numbered traces: blat DNA, -fine, -fastMap, dnax, gfServer/gfClient | `references/control-flow.md` |
| Index build, memory formula, repMatch, tile lookup, clumping, oneOff, seeding probability math | `references/index-and-seeding.md` |
| ffScore, gap penalties, minScore/minIdentity, blast e-values | `references/scoring-and-filters.md` |
| psl/pslx/blast8/blast9 fields, headers, -noHead, index file format | `references/io-formats.md` |
| Gotchas: 16-alignments-per-chromosome cap, repMatch dropout, -oneOff parse, maxIntron chaining | `references/pitfalls.md` |
| Build, tests, debugging | `references/build-run-test.md` |
| `notes`: tuning BLAT for diverged high-copy repeats (SINE use case) | `references/repeat-search-notes.md` |

Read `pitfalls.md` and `repeat-search-notes.md` first when the task is a repeat search.

## Ten facts worth knowing before touching the code

1. Defaults (DNA/DNA): tileSize 11, stepSize = tileSize, minMatch 2, maxGap 2, minScore 30,
   minIdentity 90, repMatch 1024 (`blat/blat.c:33-55`, `blat.c:714-737`, `jkOwnLib/genoFind.c:2751-2752`).
   Protein-like (`-t=dnax`/`-prot`): tileSize 5, minMatch 1, maxGap 0, minIdentity 25 (`blat.c:688-710`).
2. DNA tileSize must be 6..18; above 12 the index becomes "segmented" (12-mer head hash + 2-byte
   tail) and **repMatch filtering is disabled** (`genoFind.c:617-630`, `656-675`).
3. A tile seen >= repMatch times (saturating count) is dropped from the index entirely
   (`genoFind.c:744-769`, `841-878`, `992-1008`). Default repMatch is a table scaled by
   `tileSize/stepSize` (`genoFind.c:2716-2769`). For Alu-like SINEs nearly every consensus 11-mer
   exceeds 1024 copies in a mammalian genome, so defaults index almost none of them.
4. Query is scanned at every base with a rolling 2-bit hash; hit diagonal = tStart + qSize - qStart
   (`genoFind.c:1964-2020`). Clumps: 64 kb buckets, sort by diagonal, break when diagonal jumps
   > maxGap, keep runs with >= minMatch hits, then split if consecutive hits are > 300 bp apart
   on target (`genoFind.c:1877-1961`, `1814-1875`).
5. Alignment stage (default) = fuzzyFind `ffFind` with `ffCdna` stringency on clump +-250..500 bp
   (`jkOwnLib/gfBlatLib.c:638-667`, `414-462`); fuzzyFind self-describes as "not good for finding
   distant homologies" (`inc/fuzzyFind.h:8-11`).
6. Score = +1 match, -1 mismatch, N ignored, minus gap penalty 2 + log2(tGap)/2 + log2(qGap)
   (`lib/ffScore.c:19-72`, `lib/dnautil.c:943-969`). `-minScore` is applied in `saveAlignments`
   (`gfBlatLib.c:492-514`) and capped at 20 inside the stitcher (`jkOwnLib/supStitch.c:775-776`).
7. Identity filter: `id = 1000*(match+repMatch-2*qGapCount)/(match+repMatch+mismatch) >= minIdentity*10`,
   target gaps are free under ffCdna (`jkOwnLib/gfOut.c:131-135`).
8. **Hard cap**: standalone blat merges all bundles per target sequence name and `ssStitch` returns
   at most `ssAliCount = 16` alignments per (query, strand, target sequence)
   (`gfBlatLib.c:22`, `1330-1354`, `1555-1562`, `supStitch.c:837-841`).
9. gfServer sends at most `-maxDnaHits` (default 100) clumps per query and never supports
   `-oneOff` (`gfServer/gfServer.c:68`, `351-389`, `genoFind.c:1050-1051`).
10. `-oneOff` is parsed with `optionExists`, so `-oneOff=0` still enables it (`blat.c:717`).
    `-t=dnax` translates stop codons to 'Z', which is not a valid amino-acid tile letter, so every
    5-aa window containing a stop is unindexed/unsearchable (`lib/dnaseq.c:137-143`, `genoFind.c:728-742`).

## Open questions / unverified

- `genoFind.c:134-140` sizes segmented `endList` entries as 3 parts under `GFSERVER_HUGE` and 5
  otherwise, but `endListSetEntryValues` (`genoFind.c:151-165`) writes 4 parts for HUGE and 3
  otherwise. Looks inverted (possible overflow for blatHuge with tileSize > 12). Not run, not reported upstream here.
- Index build / query timings are from docs (`gfServer/README.blat`, `gfServer/design.txt`,
  `gfServer.c:169-172`), not measured in this session. Nothing was compiled or executed.
- Probability tables in `repeat-search-notes.md` are my own model (independent substitutions, no
  indels), not from the code or the BLAT paper.
- `src/utils` (pslReps, pslCDnaFilter, faToTwoBit, pslCheck) and `src/hg` are outside the sparse checkout.

## Change log

- 2026-10-11: initial digest at `4af84d3` (sparse checkout, read-only; focus on index/seeding/scoring
  for a diverged-repeat pre-filter use case).
{% endraw %}

## Reference pages

- [architecture](architecture/) — Architecture
- [build-run-test](build-run-test/) — Build, run, test
- [control-flow](control-flow/) — Control flow
- [data-structures](data-structures/) — Core data structures
- [index-and-seeding](index-and-seeding/) — Index, tiles, repMatch, hit lookup, clumping
- [io-formats](io-formats/) — Input and output formats
- [pitfalls](pitfalls/) — Pitfalls and gotchas
- [repeat-search-notes](repeat-search-notes/) — Notes: using BLAT as a pre-filter for diverged, high-copy repeats (SINE consensus vs genome)
- [scoring-and-filters](scoring-and-filters/) — Scoring, gap penalties, and output filters
