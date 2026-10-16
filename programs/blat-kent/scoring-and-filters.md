---
title: "blat-kent: Scoring, gap penalties, and output filters"
type: reference
generated: true
---

> Part of the [blat-kent](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/blat-kent/references/scoring-and-filters.md)

{% raw %}
# Scoring, gap penalties, and output filters

Paths relative to `src/`. Verified at `4af84d3`.

## Base-level scoring

| Function | Rule | Where |
|---|---|---|
| `dnaOrAaScoreMatch(a,b,size,match,mismatch,ignore)` | sum over positions: ignore char contributes 0 | `lib/dnautil.c:943-961` |
| `dnaScoreMatch` | +1 match, -1 mismatch, 'n' ignored | `dnautil.c:963-969`; same as `ffScoreMatch` `lib/ffScore.c:19-38` |
| `aaScoreMatch` | +2 match, -1 mismatch, 'X' ignored | `dnautil.c:979-983` |
| `dnaScore2`/`aaScore2` | per-letter versions used by HSP extension | `dnautil.c:971-977` (aa: X -> 0, match 2, else -1) |
| `bioScoreMatch` | dispatch used by the stitcher | `jkOwnLib/supStitch.c:113-124` |

There is no substitution matrix for DNA. Transitions and transversions score the same.

## Gap penalties (`lib/ffScore.c`)

`ffCalcGapPenalty(hGap, nGap, stringency)` (`:117-132`), `hGap` = target gap, `nGap` = query gap:

- `ffCdna` -> `ffCalcCdnaGapPenalty` (`:40-72`):
  `acc = 2`; if `hGap > 400000`: `+= (hGap-400000)/3000`, and if `hGap > ffIntronMax` `+= (hGap-ffIntronMax)/2000`;
  negative hGap (jump back on target) -> `hGap = -8*hGap`, squared if > 48; negative nGap -> `acc += -nGap`,
  nGap = 0; `acc += log2(hGap)/2`; if `nGap != 0` `acc += log2(nGap)` else if `hGap > 30` `acc -= 1`.
  Examples: 1 bp target insert = 2+0 = 2; 1 bp query insert = 2+1 = 3; 10/10 = 2+1+3 = 6; a 100 kb
  target-only gap = 2+8-1 = 9; a 500 kb gap = 2+33+9-1 = 43.
- `ffTight` -> `calcTightGap` (`:74-93`): `8 + (hGap-overlap) + (nGap-overlap) + overlap` (linear).
- `ffLoose` -> `calcLooseGap` (`:95-114`): `8 + ln(hGap-overlap+1) + ln(nGap-overlap+1)`.

`ffScoreSomething` (`:170-195`) = sum of block `dnaScoreMatch`/`aaScoreMatch` minus gap penalty between
consecutive blocks. `ffScore` = DNA (`:197-201`), `ffScoreCdna` = `ffScore(ali, ffCdna)` (`:203-208`,
"Perfect match is number of bases in needle"), `ffScoreProtein` (`:210-214`).

`ffScoreIntron` (`jkOwnLib/ffAliHelp.c:159-183`) is **not** part of the alignment score: it scores how close a
gap's ends are to GT/AG (CT/AC) consensus (0..4) and is used only to slide gaps (`slideIntron :186`) and in
`-fine` splice refinement (`ffSeedExtend.c:798-823`, `:1114-1131`).

## Which stringency where

| Path | Stringency | Where |
|---|---|---|
| blat dna/dna (default, -fine, -fastMap) and gfClient | `ffCdna` | `gfBlatLib.c:660-661`, `:1557`, `:588-589` |
| prot/prot | `ffTight` | `gfBlatLib.c:928` |
| -t=dnax -q=prot | `ffCdna` | `gfBlatLib.c:986` |
| -t=dnax -q=dnax | `ffLoose` (query is DNA) | `gfBlatLib.c:1292`, `:1575` |
| -t=dnax -q=rnax | `ffCdna` | same |
| fuzzyFind internal tile probabilities | per stringency `tileStrinProbMult` 0.0001/0.0005/0.0005/0.5 | `fuzzyFind.c:1345` |

## Where `-minScore` is applied

1. `ssStitch(bundle, stringency, minScore, maxToReturn)`: `if (minScore > 20) minScore = 20` (`supStitch.c:775-776`)
   so the stitcher keeps chains scoring >= 20 (or less if user set lower); stops pulling chains when one
   falls below (`:824-835`).
2. `saveAlignments` (`gfBlatLib.c:492-514`): `scoreAli` (`:464-490`, same +1/-1 and `ffCalcGapPenalty`) must
   be `>= minMatch` where the parameter is the user's `minScore` (`gfLongDnaInMem` passes `minScore`, `:1560`).
   Only then `out->out(...)` is called. Note `rangeScore` (`:931-946`) is a separate helper (gap = -4).

## `-minIdentity` (psl/pslx) (`jkOwnLib/gfOut.c:128-136`)

```
gaps = nInsertCount + (stringency == ffCdna ? 0 : hInsertCount)      # count of gaps, not bases
id   = roundingScale(1000, matchCount + repMatch - 2*gaps, matchCount + repMatch + mismatchCount)
emit iff id >= minIdentity (ppt; blat passes minIdentity*10 at blat.c:589, gfClient round(minIdentity*10))
```
`matchCount`/`repMatch`/`mismatchCount`/`countNs` are counted over aligned block bases (`:76-106`); N on either
side counts as N, not mismatch. Under `ffCdna` (all DNA/DNA modes) target gaps are free, so a spliced hit
with 100 % identity in blocks passes regardless of intron count. Default 90 % means a 300 bp block alignment
can have at most ~30 mismatches (fewer with query gaps, each costing 2 in the numerator).

For axt/maf/blast formats `minGood` is stored (`gfOut.c:444`) but `saveAxtBundle` (`:241-320`) does not test it;
the `blast` (not blast8/9) writer receives `minIdentity` (`:518-533`, `lib/blastOut.c:821-823`); whether it
filters is inside `ncbiBlastOut` (not read). Treat non-psl identity filtering as unverified.

## Blast-style scores (`lib/blastOut.c`)

Alignments are converted to `axt` with `axtScoreDnaDefault` (blastz matrix) (`gfOut.c:310-313`), then
blast8/9 print `100*matches/symCount`, `symCount`, mismatches, gap opens, 1-based coords, and
e-value = `3.0e9 * 2^-(0.0205*score)` (`blastOut.c:115-121`, `:790`), bits = `round(0.0205*score)` (`:103-107`).
The 3.0e9 database size is hard-coded; e-values are not computed from the actual database.

## Why default BLAT misses 70-80 % identity hits (verified mechanisms)

1. Seeding needs two exact 11-mers on the target grid within `maxGap` diagonal shift
   (`genoFind.c:1924-1925`, `:1938`); see probability table in `index-and-seeding.md` (~0.3 for 300 bp at 75 %).
2. Dropped overused tiles remove most seeds for high-copy repeats (`genoFind.c:857`).
3. fuzzyFind (`ffCdna`) is designed for "somewhat noisy cDNA" and gives up when the exact-tile skeleton is
   improbable (`fuzzyFind.c:1222-1229`; header `inc/fuzzyFind.h:8-11`). It extends inexactly only while
   "four of the next five bases match" (`fuzzyFind.c:34-41` comment).
4. `-minIdentity=90` rejects anything below 90 % block identity (`gfOut.c:134-135`); at 75 % the alignment is
   computed and then discarded. `-minScore=30` is rarely the binding constraint for 300 bp (225-75 = 150).
5. At most 16 alignments per (query, strand, target sequence) survive the final stitch
   (`gfBlatLib.c:22`, `:1557`, `supStitch.c:837-841`).

BLAT's own stated design targets (usage/docs, not code): mRNA/genome at >= 95 %, translated >= 80 %
("dnax ... it is best to run this on RepeatMasked data", `gfServer.c:144-145`). `notes`: the 80 % figure
follows from `-t=dnax` defaults: 5-aa tiles (15 nt), minMatch 1, so any one exactly conserved pentapeptide in
frame seeds an alignment, and protein-space scoring (+2/-1) with `-minIdentity=25` is far more permissive;
but 'Z' stops break tiles (`index-and-seeding.md`), so this trick is weak for non-coding repeats.
{% endraw %}
