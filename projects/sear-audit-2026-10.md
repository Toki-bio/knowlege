---
title: "sear audit (sear2k repo and its callers), with an L1 test of minimap2 mode"
type: notes
---

> Type: notes. Every finding below marked **verified** was reproduced in a
> container on 2026-10-12; the commands are in this session's scratch work and
> can be rerun with `projects/sine-search-bench/`. Code references are to
> Toki-bio/sear2k `sear` at `85421e7` unless the vendored SINEderella copy
> (`90d8ac4`) is named.

## Status (2026-10-12)

**minimap2 mode fixed** on `main` of Toki-bio/sear2k (`923ffb6`) and Toki-bio/SINEderella
(`38d811a`), and merged into both `fast-mode` branches: query indexed, genome streamed
as overlapping windows (`seqkit sliding`, step max(len, 2 kb), width step + len + 1 kb)
into `minimap2 -x map-ont -k 11 -w 3 -c --secondary=no`; identity = PAF `$10/$11`,
score = `AS:i`; cascade removed; output prefix from the file name; `sear` executable in
sear2k. ssearch36 mode now warns and skips when len x Slf > 2500. The other ssearch36
findings (1, 2, 3, 5-8) are **not** fixed.

Verified after the fix (chr21, RepeatMasker truth, 4 threads):

| query | old sear | fixed sear | blastn dc-megablast, same filters |
|---|---|---|---|
| L1.3 (6,059 bp), full-length L1 kept by the 0.8 length rule (of 73) | 4 | 52 (65 at any length) | 50 |
| L1PA7 genomic copy (6,468 bp), loci | – | 62 | 47 |
| L1.3 3' end (1,500 bp), loci | 1 | 305; 89 % of dc-megablast's 308 | 308 |
| identity reported for 98 % copies | ~70 % | 92-98 % | |
| runtime, 6 kb query | 2 s (chr21) | 9 s (chr21), 106 s (1.65 Gbp lizard genome) | 4 s (chr21) |

Every locus the fixed mode reports lies in an annotated L1. For reference, nhmmer
finds 72 of 73 full-length L1s but takes 21 min on chr21 for one 6 kb query, and
ssearch36 on 8 kb windows with the 0.8 length rule keeps 58 in 7 min: the fixed
minimap2 mode is within that range at 9 s. Seeds `-k 11 -w 3` beat
map-ont's default `-k 15 -w 10` (1.5 kb query: 89 % vs 77 % of dc-megablast's loci) and
`asm20` (35 of 73). SINE searches (ssearch36 path) reproduce all loci of a pre-fix run.
The 23 full-length L1s an L1HS query does not keep are other subfamilies (L1PB, L1MA,
L1PA10-16) whose alignable span is under 0.8 x 6 kb; dc-megablast drops the same ones.

## Summary (as audited, before the fix)

- **minimap2 mode (queries ≥ 1000 bp) does not work for repeats.** On human
  chr21 with the L1.3 element (6,059 bp) it found **4 of 73** full-length L1s
  that RepeatMasker annotates; blastn found 68. Its identity column reports
  ~70 % for copies that are 98 % identical, and its bitscore column is always 0.
- **ssearch36 mode is sound for SINE-length queries**, but has five traps that
  silently lose or corrupt results under conditions SINEderella avoids and
  manual use does not.
- The **sear2k repo copy is behind the SINEderella copy**; two older scripts in
  the repo have different option semantics.

## minimap2 mode, tested with L1 on chr21

Setup: GRCh38 chr21 (46.7 Mbp, Ensembl), RepeatMasker annotation from the UCSC
API (12,508 L1 rows; 73 full-length ≥ 4.8 kb; 307 young ≥ 1.8 kb with
divergence ≤ 15 %). Query L1.3 (GenBank L19088, an L1HS). A copy counts as
found when hits cover ≥ 50 % of it. 4 threads.

| method | time | full-length L1 found (73) | young ≥ 1.8 kb found (307) |
|---|---|---|---|
| `sear` (minimap2 mode, Slf 0.8 or 0.3) | 2 s | **4** | 4 |
| `sear --force-ssearch` | 176 s | **0** | 0 |
| minimap2 `-x asm20 -c -f 0 -N 1e5 -p 0.1`, genome indexed | 1.6 s | 10 | 17 |
| minimap2 `-x map-ont -c -f 0`, genome indexed | 2.4 s | 11 | 20 |
| **minimap2 roles swapped**: L1 indexed, genome as 10 kb windows (7 kb overlap), `-x map-ont -c` | 5.7 s | **65** | **264** |
| blastn `-task megablast` | 0.8 s | 56 | 205 |
| blastn `-task dc-megablast` | 3.6 s | 68 | 284 |
| blastn `-task blastn` | 3.8 s | 68 | 288 |
| nhmmer, L1.3 as query | 1,258 s | 72 | – |
| ssearch36, 8 kb windows / 2 kb overlap, 0.8 length rule (what `--force-ssearch` would need) | 421 s | 58 | – |

The four loci `sear` reports are all real L1PA2 copies; nothing it reports is
outside an L1. It just misses almost everything.

Why (from the minimap2 digest and the tests):

1. **Seeds of a high-copy repeat are discarded.** With the genome indexed,
   nearly every minimizer of an L1 query occurs hundreds of times and falls
   above minimap2's occurrence cutoff (`-f`); per streak of repetitive seeds
   only about one is rescued (`seed.c:56-96`). `-f 0` and smaller k/w did not
   help (10-11 of 73). Swapping roles makes the index one sequence, so no
   filter can fire.
2. **No `-c` (`sear:323-326`)**: without base-level alignment, PAF column 10 is
   a seed-based count, so `pct_id = $10 / min(qspan, tspan)` (`sear:337-340`)
   reports 68-73 % for copies minimap2 `-c` measures at 97.9-98.2 %
   (RepeatMasker 1.2-1.8 % divergence). With the default Homology 65, copies
   that really are 80-85 % identical would be filtered out even when found.
3. **`AS:i` exists only with `-c`**, so column 7 is 0 for every row
   (`sear:342-345`). `-b` (best N by column 7) and `-s` (top N) then pick
   arbitrary loci.
4. **The cascade never escalates.** It stops at the first preset with *any* raw
   alignment, before filtering (`sear:349-354`), so asm20 runs only if asm10
   finds nothing at all. The comment and README say asm5 → asm10 → asm20; the
   code runs asm10 → asm20 (the SINEderella copy fixed the comment, not the
   logic).
5. **`--force-ssearch` cannot work for long queries**: windows are 2.5 kb, so
   for a query longer than 3,125 bp at Slf 0.8 no hit can pass the length rule.
   It scans the whole genome and returns nothing, without a warning.

What would fix it (proposal, not applied): in minimap2 mode, index the query
and stream the genome as overlapping windows at least query-length + a few kb
wide, run with `-c`, take identity as `1 - NM/alignment length` and the score
from `AS:i`, keep `-p`/`-N` permissive, and merge window hits back with the
existing window-offset code. Or call blastn (`-task dc-megablast` or `blastn`)
on the same windows or with a database. Both gave ~90 % recall in 4-6 s on
chr21. Make `--force-ssearch` refuse queries longer than 2,500 × Slf.

## ssearch36 mode findings

| # | finding | severity | verified | who is exposed |
|---|---|---|---|---|
| 1 | **Splits of another genome are reused.** The check is `compgen -G "*.2k.part_*.bnk"` (`sear:397`), any genome. With `-k`, or leftovers, a second genome in the same directory is never searched: its run scanned the first genome's parts, wrote nothing under its own prefix, and reported 0 loci. | high | yes | manual use; SINEderella and the walkers use one genome per directory |
| 2 | **Cleanup globs delete unrelated files in the working directory**: `minus_*`, `plus_to_add_*`, `[0-9]*_*.bed`, `*_s.fai` at start and end (`sear:268-288`), and without `-k` every `*.2k.bnk*` / `*.2k.part_*` (`sear:683`). A file `2_mydata.bed` in the run directory was deleted. | high | yes | manual use |
| 3 | **Long contig names corrupt coordinates.** The genome offset is parsed from the window name `contig_sliding:start-end` as ssearch36 prints it in `-m 8C` column 2. That column is cut at the `-w` width: 60 characters by default, 95 with `sear`'s `-w 95`. Names whose window ID exceeds ~95 characters lose the offset; with 82-character contig names (after `_` → `@U@`, +2 per underscore) the run gave 30 rows with impossible spans instead of 408. SPAdes `NODE_1_length_3000000_cov_25.123456` stays under the limit. | medium | yes | assemblies with very long headers |
| 4 | Query longer than 3,125 bp in ssearch36 mode returns nothing (see above). | medium | yes | `--force-ssearch`, or `minimap2` missing from PATH (sear then falls back to ssearch36 silently, `sear:298`) |
| 5 | `-b` best-hit loci are shifted 1 bp right relative to the main BED: the main path applies `bedtools slop -l 1 -r -1` (`sear:549`) to undo the 1-based window start, the best-hit path does not (`sear:564-582`). 5 of 5 checked. | low | yes | `-b` users |
| 6 | A multi-sequence query fails with `awk: cannot open "BEGIN{...}"`, exit 2, no output (`sear:273-274` computes one length per record). | low | yes | anyone passing a bank as query (that is what `sear_multi` is for) |
| 7 | `-s N` keeps the first N hits in part order, not the best N. LINE_walker and SINE_walker call `sear -s 50` and then take their "top 20 by bitscore" from those 50. | design | by reading | the walkers |
| 8 | Intermediates have fixed names in the working directory; two `sear` runs in one directory overwrite each other. | low | by reading | parallel manual runs |
| 9 | ssearch36 flags: `-g -3` is gap-extension −3, not "both strands" (both strands are the default anyway); a single `-E 2` makes the secondary-HSP threshold 1.0; no `-Z`, so E-values depend on part size. | low | earlier session | everyone |
| 10 | (sear2k repo only) In minimap2 mode the output prefix is the first 3 characters of the genome *path* (`sear:279`): with `/tmp/...` it wrote `/tm-q.bed` to the filesystem root. The SINEderella copy uses `basename`. | medium | yes | sear2k users with absolute paths |

## Repository state

- `sear2k/sear` and `SINEderella/sear` have diverged in 14 hunks. The
  SINEderella copy carries fixes the upstream lacks (basename prefix, raw-hit
  file for the satellite stage, cascade comment). `SINEderella-dev/sear` is
  older still. One copy should be the source and the others should vendor it.
- `sear2k` and `sear2k_bak` are older scripts: `-T 24` hard-coded,
  `-m` *skips* minus rounds (the opposite of `sear`), default Slf 0.9. The
  README documents only `sear`. Archive them or mark them deprecated.
- `annotate_aligned_by_ssearch_seqkit.sh`: `xargs -n1 -P N -I {}` mixes
  exclusive options (GNU xargs lets `-I` win, so it works) and builds a shell
  string from the database path, which breaks on paths with spaces or quotes.
  Minor.

## Callers

| caller | how it calls `sear` | exposure |
|---|---|---|
| SINEderella step 1 and `--add` | `sear -k query genome.clean.fa 0.8 65 50`, one consensus at a time, own `searches/` directory per run, splits deleted after the last query | ssearch36 path, safe from 1, 2, 8. **Any consensus ≥ 1000 bp goes silently to minimap2 mode** (finding table above), e.g. a composite or LINE-derived entry in a bank. |
| SINEderella `--fast` (branch `fast-mode`) | same, `SEAR_MODE=fast` | same exposure; the patch does not touch minimap2 mode |
| LINE_walker / SINE_walker | `THREADS=n sear -s 50 -k query genome.fa` in a shared work directory per genome, parse `gen-<query>.bed` sorted by bitscore | finding 7 (first-N, not best-N). If a walk step's query reaches 1000 bp, it switches to minimap2 mode, where bitscore is 0 and the "top hits" are arbitrary. |
| Tal project | through SINEderella | as SINEderella |

## Corrections to earlier notes

- `-w 95` is not display-only in `-m 8C` output: it sets the width of the
  subject-ID column (60 by default). It is what keeps `sear`'s window offsets
  intact for normal contig names. Any wrapper that drops `-w` (my benchmark
  scripts did, harmlessly, because their IDs were short) must keep IDs under
  60 characters.
