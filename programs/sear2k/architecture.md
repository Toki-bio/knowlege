---
title: "sear2k: Architecture"
type: reference
generated: true
---

> Part of the [sear2k](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/sear2k/references/architecture.md)

{% raw %}
# Architecture

Full commit SHA: `85421e7964178974fe856a443508560613f88aca`
(2026-03-05, "CRITICAL FIX: filter on merged span (\-\) not max fragment length
- was checking wrong column"; the clone is shallow, so this is the only commit
visible).

Repo URL: https://github.com/Toki-bio/sear2k

## Directory map

| Path | Lines | Role |
|---|---|---|
| `sear` | 653 | Current tool. `set -euo pipefail`, option parser, ssearch36 fragment mode + minimap2 mode, `-b` MAFFT mode, `--slop`, `-k`, `--force-ssearch`. This is what SINEderella's `step1_search_extract.sh` calls (as a vendored copy). |
| `sear2k` | 370 | Previous generation. Same ssearch36 call with hard-coded `-T 24`; minus-bank rounds ON by default (`-m` turns them off); Slf default 0.9; no minimap2, no `--slop`. |
| `sear2k_bak` | 200 | Original v0.2 (`sear2k_bak:4`). Positional args only, no bitscore column, same split and ssearch36 call. Kept as history. |
| `annotate_aligned_by_ssearch_seqkit.sh` | 109 | Unrelated helper: for every sequence of an alignment, `ssearch36 -m 8` against a consensus DB (one process per sequence via `xargs -P`), then prefixes `topHit=...|bitscore=...` to the aligned headers. |
| `README.md` | 147 | Describes `sear` (not `sear2k`): options, minimap2 cascade, BED7 format. |

No tests, no data, no Makefile, no CI.

## Module dependency sketch

```
sear (bash)
 ├── samtools faidx           genome index; also per-part .fai for minus rounds
 ├── seqkit sliding/split2    2500 bp windows, step 2000, 20 parts         (ssearch36 mode)
 ├── ssearch36 -m 8C          query vs each part                           (ssearch36 mode)
 ├── minimap2 -x asm10/asm20  query vs whole genome                        (query >= 1000 bp)
 ├── bedtools merge/sort/slop/complement/getfasta
 └── mafft --localpair        only with -b/--best
```

Everything is a pipeline of awk filters between these tools; there are no
functions except `usage`, `cleanup_query` (`sear:254-260`) and, in the
SINEderella copy, `to_genomic`.

## Entry points

| Command / mode | Entry | Notes |
|---|---|---|
| `sear [opts] query.fa bank.fa [Slf=0.8] [Homology=65] [Flank=50] [KeepSplits=0]` | `sear:62-139` option loop, `sear:216-219` positionals | mode chosen at `sear:265-272` |
| ssearch36 fragment mode | `sear:354-499` | default for query < 1000 bp or `--force-ssearch` |
| minimap2 mode | `sear:274-352` | query >= 1000 bp and `minimap2` on PATH |
| best-hit + MAFFT | `sear:503-575` | `-b [N]` |
| final extraction (`.bed` -> `.bnk`) | `sear:577-583` | always |
| `sear2k [opts] query bank [Slf=0.9] [Homology=65] [Flank=50] [KeepSplits=0]` | `sear2k:37-83`, `sear2k:146-149` | search loop `sear2k:175-241` |

## Version lineage (what changed, where)

| Feature | `sear2k_bak` | `sear2k` | `sear` (this repo) | `sear` in SINEderella `90d8ac4` |
|---|---|---|---|---|
| ssearch36 flags | `-T 24` (`:58`) | `-T 24` (`:190`) | `-T ${THREADS:-nproc}` (`:396`) | same (`:404`) |
| Slf default | 0.9 (usage `:10`) | 0.9 (`:146`) | 0.8 (`:216`) | 0.8 (`:216`) |
| bitscore column in BED | no (6 cols, `:99`) | yes, col 7 (`:258`) | yes | yes |
| minus-bank rounds | always | default on, `-m` skips (`:54-56`) | default off, `-m` enables (`:81-83`) | same |
| `-s/--stop`, `-b/--best` | no | yes | yes | yes |
| `-k`, `--slop`, `--force-ssearch` | no | no | yes | yes |
| minimap2 mode | no | no | yes; text says asm5->asm10->asm20, loop runs asm10, asm20 (`:36,293`) | text fixed to asm10 -> asm20 (`:36-37,296`) |
| output prefix | `cut -c1-3` of bank path | same | same (`:251`; bug for absolute paths in minimap2 mode) | `basename` first (`:254`) |
| raw-hit file (`*.rawhits.tsv`, every merged hit at the identity cut, any length) | no | no | no | yes (`:341,387,417,540-544`) |
| empty-ssearch guard before `bedtools merge` | no | no | yes (`:404-411`) | yes |
| `_sliding` truncated-suffix handling (`_slidin`, `_slidi`, ...) | no | no | yes (`:463-472`) | yes |

The SINEderella copy is the live one; this repo's `sear` is one step behind
(2026-03 vs 2026-10). Fixes flow from SINEderella back here only by hand.

## Design decisions the authors made (and why, if stated)

- **Fragment the genome, not the query.** Comments: "splitting the genome into
  overlapping fragments to avoid boundary misses" (`README.md:3`). Windows of
  2500 bp with a 500 bp overlap mean a <= 500 bp SINE copy is whole in at least
  one window. `notes`: this also keeps each library sequence short, which is
  what makes ssearch36's per-sequence Smith-Waterman tractable and what makes
  `-E 2` meaningful per window.
- **Merge before filtering.** Several local alignments of one copy (gapped,
  split by indels) are merged within the window and the span of the merged
  interval is what must reach 0.8 x query length (`sear:405-407`). The last
  commit says this check previously read the longest fragment, not the span.
- **minimap2 only for long queries.** "asm5 is too strict for diverged copies"
  (SINEderella `sear:36-37`); below 1000 bp minimap2 is not used at all.
  SINEderella's `docs/MASKING.md:34` records that minimap2 `-k 11 -w 5`
  returned nothing for SINE-sized families (repeat k-mers exceed `mid_occ`).
- **Header rewriting instead of escaping.** Underscores are replaced by `@U@`
  (shorter than the earlier `@UNDERSCORE@`, comment `sear:194`) because the
  window-offset parser splits on `_` and ssearch36 truncates long identifiers.
{% endraw %}
