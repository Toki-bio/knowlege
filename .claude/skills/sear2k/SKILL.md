---
name: sear2k
description: Internals of sear2k / sear (Toki-bio), the user's own Bash tool that finds short sequences (SINE consensuses, 100-500 bp) in a genome by splitting the genome into 2.5 kb sliding windows and running ssearch36 on each 1/20th part, then merging hits back to genome BED and extracting copies with flanks. Use when reading, modifying, benchmarking or replacing the ssearch36 search step (its exact ssearch36 / minimap2 command lines, hit parsing, filters, minus-bank rounds, best-hit MAFFT mode) or when deciding how to accelerate it.
---

# sear2k / sear internals

> Type: reference (facts verified against the pinned commit) with marked `notes`.

## Pinned source

| Field | Value |
|---|---|
| Repo | https://github.com/Toki-bio/sear2k |
| Version / tag | none; `sear2k` header says "Version: 0.2 + bitscore + -b/--best ..." (`sear2k:4`) |
| Commit | `85421e7` (full SHA in references/architecture.md) |
| Digested on | 2026-10-11 |
| Language / build | Bash (3 scripts, 1 backup), no build; calls ssearch36, seqkit, bedtools, samtools, mafft, minimap2 |
| License | MIT (`README.md:146`, `sear2k:5`) |

## What it does in three sentences

`sear` takes one query FASTA (a SINE consensus) and a genome FASTA, cuts the
genome into 2500 bp windows stepping 2000 bp (`seqkit sliding`) and splits
those windows into 20 part files, then runs `ssearch36 -m 8C` (full
Smith-Waterman, both strands, shuffle statistics `-z 11`, `-E 2`) of the query
against each part in turn. Hits are merged per window with `bedtools merge -s`,
kept when identity >= 65 % and merged span >= 0.8 x query length, mapped back
to genome coordinates by parsing the window offset out of the `_sliding:`
header, merged again, and written as BED7 plus a flanked FASTA (`.bnk`). For
queries >= 1000 bp it switches to `minimap2 -x asm10/asm20` on the whole
genome; `-b` keeps the top-N hits by bitscore and aligns them with MAFFT.

## Where to look (read on demand)

| Question | File |
|---|---|
| Which script is which, version lineage, who calls whom | `references/architecture.md` |
| Exact ssearch36 / minimap2 / mafft command lines, flag semantics, how `-m 8C` output is parsed (Q1) | `references/aligner-invocation.md` |
| Step-by-step trace of ssearch36 mode, minus-bank rounds, minimap2 mode, best mode (Q2, Q3) | `references/control-flow.md` |
| Inputs, BED7/bnk outputs, intermediates, naming | `references/io-formats.md` |
| Where the wall-clock goes, parallelism, author's speed notes (Q4) | `references/performance.md` |
| Gotchas: `-g -3` parsing, taxname path bug, cascade text mismatch, `-E` relative to part size | `references/pitfalls.md` |
| Run it, and what tests exist (none in this repo) (Q5) | `references/build-run-test.md` |

Read only the file that answers the current question. For the strategy question
"replace or accelerate ssearch36" read `aligner-invocation.md` then
`performance.md`; the sibling skill `sinederella` holds the measured
blastn/minimap2/ssearch comparison and the real-genome timings.

## Ten facts worth knowing before touching the code

1. Three generations live side by side: `sear2k_bak` (v0.2, 200 lines),
   `sear2k` (v0.2 + flags, 370 lines), `sear` (653 lines, `set -euo pipefail`,
   minimap2 mode). SINEderella vendors a newer `sear` (702 lines) with a
   raw-hit file; see `references/architecture.md`.
2. The ssearch36 call is identical in all three:
   `ssearch36 -g -3 -T $T -Q -n -z 11 -E 2 -w 95 -W 70 -m 8C query part 3`
   (`sear:396`, `sear2k:190`, `sear2k_bak:58`). Only `-T` changed from a
   hard-coded 24 to `${THREADS:-nproc}`.
3. Library = one of 20 part files of 2.5 kb windows, never the whole genome
   (`sear:372-373`). Query = one consensus. The 20 parts are searched
   sequentially; parallelism is only ssearch36's own `-T`.
4. Filters are applied after `bedtools merge -s` within a window: identity
   (col 7 after merge, `max` of hit identities) >= Homology (65) and merged
   span >= int(Slf x querylen) with Slf = 0.8 (`sear:405-407`; 0.9 in `sear2k:198`).
   The last commit's message says this span check used to read the wrong column.
5. Window -> genome mapping parses the `_sliding:start-end` suffix, sums the
   positive numbers, adds to the hit coordinates, then `bedtools slop -l 1 -r -1`
   and a second strand-aware merge (`sear:456-499`). Two overlapping windows
   therefore collapse a boundary-straddling copy into one locus.
6. `-m/--minus` enables repeated rounds per part: hits are complemented out of
   the part with `bedtools complement | getfasta -s` and the remainder is
   searched again until a round returns 0 hits (`sear:434-449`). In `sear2k`
   the default was the opposite: rounds on, `-m` to skip (`sear2k:54-56,221-223`).
7. minimap2 mode (`sear:266-352`): query >= 1000 bp and minimap2 on PATH;
   `minimap2 -x asm10|asm20 --secondary=yes -p 0.3 -N 10000 genome query`,
   PAF parsed in awk, identity = matches / min(qspan, tspan); stops at the
   first preset with any raw alignment. The usage text still says asm5 (`sear:36,275,291`).
8. `-s/--stop N` (default 1000) stops scanning further parts once N filtered
   hits have accumulated (`sear:419-424`); in minimap2 mode it keeps the top N
   by AS instead (`sear:341-347`).
9. `-b/--best [N]` re-ranks merged genome hits by bitscore (col 7), takes N,
   slops 50 bp (or `--slop`), `getfasta -s -nameOnly`, and runs
   `mafft --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder` with query
   prepended (`sear:503-575`).
10. Genome headers with `_` or tabs are rewritten (`_` -> `@U@`) into a
    `*_temp_clean.*` copy because the offset parser splits on `_`
    (`sear:180-214`, restored at `sear:585-612`).

## Open questions / unverified

- FASTA36 option semantics for `-g -3` (gap-extension penalty -3 vs
  "forward strand only"), `-w 95 -W 70` (display widths, no effect with `-m 8C`),
  the trailing `3` (ktup, ignored by ssearch) are my reading of the FASTA36
  manual; ssearch36 is not installed here, so not run-verified. SINEderella's
  MANUAL reads the same flags differently (see `references/pitfalls.md`).
- No timing data exists in this repo; all wall-clock numbers come from the
  SINEderella repo (sibling skill).
- `sear2k_bak` has a `Secs_HMS` timer defined (`sear2k_bak:20-21`) but I did
  not trace whether it is printed anywhere.

## Change log

- 2026-10-11: initial digest at `85421e7`.
