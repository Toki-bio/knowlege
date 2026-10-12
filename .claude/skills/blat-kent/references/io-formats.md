# Input and output formats

Paths relative to `src/`. Verified at `4af84d3`.

## Inputs accepted

| Format | Reader | Notes / quirks |
|---|---|---|
| FASTA (target or query), detected by first char `>` | `gfClientFileArray` `jkOwnLib/gfClientLib.c:16-52`; target via `faReadAllMixed` `:188`; queries via `faMixedSpeedReadNext` `blat/blat.c:390` | A non-FASTA argument is treated as a **list file** of file names (`:48-51`). `.gz/.Z/.bz2` names are single files. Target lower case = unmasked after load (`:163-168`). Duplicate target names abort (`genoFind.c:1429-1442`). |
| `.nib` | `nibLoadAllMasked`, `blat.c:341-353`; `gfCountTilesInNib` `genoFind.c:771-794` | one sequence per file; `file.nib:start-end` subrange (`blat.c:73-79`) |
| `.2bit` (optionally `file.2bit:seq:start-end`) | `twoBitIsSpec`, `blat.c:354-385`; `gfCountTilesInTwoBit` `genoFind.c:817-839` | gfServer requires `.2bit`/`.nib` (`genoFind.c:1072-1073`); `gfServer index` needs exactly one `.2bit` (`gfServer.c:1157-1158`) with name `<base>.{untrans,trans}.gfidx` (`:1133-1152`) |
| `.ooc` | `oocMaskCounts` `jkOwnLib/ooc.c:9-53` | 32-bit signature, tileSize, tile values; tileSize must match (`:30-32`); byte-swapped files handled |
| RepeatMasker `.out` (`-mask=out`/`file.out`, `-repeats`) | `maskFromOut` `gfClientLib.c:63-...` | repeats with `percDiv+percDel+percInc <= minRepDivergence` (15) are masked (`:94-95`) |
| `.gfidx` | `genoFindIndexLoad` `genoFind.c:560-583` | mmap, checks magic/version/address size/isTrans (`:493-511`) |

Query types `-q=dna|rna|prot|dnax|rnax`, target `-t=dna|prot|dnax` (`blat.c:81-93`); `rna/rnax` converts u->t
(`blat.c:320-321`) and trims polyA by default (`:684-687`).

## Outputs produced (`-out=`, `gfOutputAny` `jkOwnLib/gfOut.c:535-577`)

| Format | Writer | Fields | Flags |
|---|---|---|---|
| `psl` (default) | `savePslx` `gfOut.c:36-194` | 21 tab columns below | `-noHead` suppresses 5-line header (`gfOut.c:469-470`, header `lib/psl.c:568-573`, `:541-549`: "psLayout version 3", blank, two label lines, dashes) |
| `pslx` | same with `saveSeq` | 21 columns + `qSeqs` + `tSeqs` (comma lists of block sequences, `:170-184`) | `-noHead` |
| `axt` | `saveAxtBundle`/`axtQueryOut` `:241-335` | blocks split at gaps > 8 (`ffNextBreak(sAli, 8, ...)` `:262`) | |
| `maf` | `mafQueryOut` `:336-358` | | |
| `sim4` | `sim4QueryOut` `:393-425` | DNA only (`:511-512`) | |
| `blast`, `wublast`, `xml` | `lib/blastOut.c` via `axtBlastOut` `:799-832` | pairwise text | uses `databaseName/SeqCount/Letters` from blat (`blat.c:585-590`) |
| `blast8` | `tabBlastOut(withComment=FALSE)` `blastOut.c:731-797` | 12 cols: qName, tName, %id (`100*matches/symCount`), aln length (symCount incl. gaps), mismatches, gap opens, qStart, qEnd (1-based), sStart, sEnd, e-value, bit score | one row per gap-free-ish axt piece (split at gaps > 8), not per psl line |
| `blast9` | same with `#` comment lines (`:740-757`) | | |

### psl columns as written (`gfOut.c:152-169`)

1 `matches` 2 `misMatches` 3 `repMatches` (only with `-mask/-repeats` maskHash, `:92-99`) 4 `nCount`
5 `qNumInsert` 6 `qBaseInsert` 7 `tNumInsert` 8 `tBaseInsert` 9 `strand` (`+`/`-` query; translated modes
append target strand, `reportTargetStrand`, `:155-156`, `blat.c:433`, `:442`) 10 `qName` 11 `qSize`
12 `qStart` 13 `qEnd` 14 `tName` 15 `tSize` 16 `tStart` 17 `tEnd` 18 `blockCount` 19 `blockSizes,`
20 `qStarts,` 21 `tStarts,`. Coordinates 0-based half-open; for `-` strand, `qStart/qEnd` are flipped to
forward-strand coordinates (`:137-144`) while `qStarts` stay in reverse-strand coordinates (`:165-166`,
standard psl convention). Translated target coords are multiplied back to nucleotides via `trans3GenoPos`
(`:72-73`, `:169`).

### Getting all hits

- There is no "report all" switch. Limits: `-minScore`, `-minIdentity`, the 16-per-target-sequence stitcher cap
  (`gfBlatLib.c:22`, `supStitch.c:837-841`), gfServer `-maxDnaHits=100` clumps per query (`gfServer.c:68`,
  `:383-384`). See `pitfalls.md` and `repeat-search-notes.md`.
- `-maxIntron` (default 750000, `inc/fuzzyFind.h:63`) is **not** irrelevant for repeat searches: with the default
  the stitcher may chain two partial copies up to 750 kb apart into one multi-block psl line (gap penalty
  only ~9-43, `lib/ffScore.c:40-72`); `cutAtBigIntrons` runs only when `-maxIntron` differs from the default
  (`supStitch.c:815-820`). Set it small (e.g. `-maxIntron=50`) to force one psl line per copy. Changelog notes
  a fixed segfault for `maxIntron <= 11` with `tileSize <= 11` (`blat/version.doc` 34x11).
- `-noHead` only drops the header; `-dots=N` prints progress to stdout.
- Post-filters mentioned in usage (`pslReps`, `pslCDnaFilter`, `pslUniq`, `blat.c:157-158`) live in `src/utils`
  and `src/hg`, which are not in this sparse checkout.

## On-disk index / intermediate files

- `.gfidx` layout: `genoFindIndexFileHdr` (`genoFind.c:462-476`) then one (`untrans`) or six (`trans`)
  `genoFindFileHdr` sections (`:196-221`) with 8-byte-aligned (`mustSeekAligned :183-193`) arrays:
  sources (file names + start/end, `:262-307`), `listSizes`, `lists` (offsets relative to block start,
  `:324-355`) or `endLists`. Written by `genoFindIndexWrite` (`:534-558`), loaded with `mmap` +
  `madvise(MADV_RANDOM|MADV_WILLNEED)` (`:560-583`). 32-bit and 64-bit (`Huge`) files are incompatible (`:503-508`).
- gfServer wire protocol: each command starts with `gfSignature()` `"0ddf270562684f29"` (`genoFind.c:64-70`);
  `query <size>` then raw sequence; server answers clump lines `qStart qEnd fileName tStart tEnd hitCount`
  and `end` (`gfServer.c:370-373`, `:917`); translated queries also send a `tileSize N` line and per-clump hit
  lists (`:404`, `:420-428`). Dynamic server naming: `$rootdir/$genomeDataDir/$genome.{2bit,untrans.gfidx,trans.gfidx}`
  (`gfServer.c:120-125`).
