# Control flow

Line numbers are for `sear` at `85421e7`.

## `sear [opts] query.fa bank.fa [Slf] [Homology] [Flank] [KeepSplits]`

### Setup

1. Parse options into `BEST, BEST_N, DO_MINUS, STOP_EARLY, STOP_N, KEEP_FLAG, FORCE_SSEARCH, SLOP_L/R` (`:51-139`); positionals `query_in bank_in` (`:146-147`), defaults `Slf=0.8 Homology=65 FlankSize=50 KeepSplits=0` (`:216-219`).
2. Require `bedtools samtools` (+ `mafft` with `-b`) (`:158-167`); ssearch36/seqkit are checked only in ssearch36 mode (`:360-365`).
3. `samtools faidx` the bank if no `.fai` (`:173-178`).
4. If any header has `_` or a tab: write `<bank>_temp_clean.<ext>` with `_` -> `@U@` in the first header word only, index it, use it as `WORKING_BANK`; reuse an existing clean copy if newer than the original (`:180-214`).
5. Delete stale intermediates (`:238-243`); register `cleanup_query` on EXIT (`:254-260`).
6. Query length `len` = total bases over all records of the query file (`:245`, so a multi-record query is summed, not per record); `flen = int(Slf*len)` (`:246`).
7. `qname` = query path without extension (`:250`); `taxname` = first 3 chars of the bank path (`:251`, overwritten per part at `:387` in ssearch36 mode).
8. Mode: `len >= 1000 && !FORCE_SSEARCH && minimap2 on PATH` -> minimap2 (`:265-272`).

### ssearch36 mode (`:354-499`)

9. If no `*.2k.part_*.bnk` in cwd: `seqkit sliding -g -s 2000 -W 2500 bank > bank.2k.bnk; seqkit split2 -p 20 -O . bank.2k.bnk` (`:369-374`). Splits are reused across queries when present (this is how `step1_search_extract.sh` amortises them with `-k`).
10. For each part (`:381`), `round=1`, `Pl=1`; while `Pl != 0` (`:392`):
    1. Run ssearch36 (`:396`), awk to window-BED (`:397-401`), `sort -V` to a temp file.
    2. If non-empty: `bedtools merge -s -c 4,5,6,7,8 -o max,max,distinct,max,max`, keep `pctid >= Homology && hitlen >= flen`, `sort -V` -> `plus_to_add_<round>_<part>.bed` (`:404-411`).
    3. `Pl` = line count; append to `bank.all_hits.txt`; `HIT_COUNT += Pl`; if `STOP_EARLY && HIT_COUNT >= STOP_N` -> `STOP_NOW=1; break` (`:414-424`).
    4. Index the part, sort its `.fai` (`:426-427`).
    5. If `Pl == 0` -> break. If `!DO_MINUS` -> break (so by default exactly one round per part) (`:429-436`).
    6. Minus bank: concatenate all `plus_to_add_*_<part>.bed`, `bedtools complement -g part_s.fai | sort | getfasta -s`, clean `:()` from headers -> `minus_<round>_<part>`; `bank=minus_...; round++` and loop (`:438-449`). Each round therefore searches the part minus everything found so far, until a round finds nothing.
11. After all parts (or STOP): compute per-hit window offsets from the header (`:456-460`): replace `_sliding:`/`_sliding_` by a space, print field 2 (`2001-4500`), turn `-` into ` -` so the end becomes the negative field `-4500`, and sum only the positive fields: the offset is the 1-based window start (2001). Adding it to the 0-based in-window start and the 1-based in-window end gives coordinates one too high, which `bedtools slop -l 1 -r -1` (`:497`) shifts back by one. `notes`: read, not run; a one-line awk test with a hand-written header is cheap insurance before reusing this parser.
12. Build `bank.stripped.genome` (chrom lengths with every `_sliding*` suffix variant stripped, max per chrom) (`:463-480`).
13. `paste all_hits all` -> add offset to start/end -> strip suffix from chrom -> print `chrom gstart gend pctid qspan strand bits` -> `bedtools sort | slop -l 1 -r -1 | merge -s -c 4,5,6,7 -o max,max,distinct,max` -> `<taxname>-<qname>.bed` (`:482-499`).

### minimap2 mode (`:274-352`)

9'. Auto `-b` size if `-b` given without N: 10 for `len < 5000`, else 5 (`:278-285`).
10'. `stripped.genome` straight from the `.fai` (`:288-289`).
11'. For `preset in asm10 asm20`: run minimap2, parse PAF to `tname tstart tend pct_id tlen strand AS` (`:293-319`); break on first non-empty output (`:321-326`).
12'. No hits -> empty BED (`:329-331`). Else `bedtools sort | merge -s -c 4,5,6,7 -o max,max,distinct,max | awk $4>=hom && ($3-$2)>=flen` (`:334-338`); with `-s`, keep top `STOP_N` by col 7 (`:341-347`).

### Best mode (`:503-575`, both modes)

13. minimap2: `sort -k7,7nr bed | head -N` (`:508-509`). ssearch36: recompute genome coordinates from `all_hits`, `bedtools sort | merge -s`, `sort -k7,7nr | head -N` (`:511-531`).
14. To BED6 with name `tmp` (`:534-536`), `bedtools slop -s -l BestSlopL -r BestSlopR` (default 50/50) (`:538-539`), name = `chrom:start-end` (`:541-543`), dedupe on `(chrom,start,end,strand)` (`:545-547`), `getfasta -s -nameOnly` (`:549-550`).
15. `cat query; echo; cat hits` -> MAFFT (`:556-564`); fail if output empty (`:566-569`); keep only `<taxname>-<qname>.best<N>.mafft.fa` (`:614-623`).

### Final extraction and cleanup (`:577-653`)

16. `bedtools slop -s -l MainSlopL -r MainSlopR` on the BED (`:581`), `getfasta -s` -> `<taxname>-<qname>.bnk` (`:582`). Flank = positional `Flank` (50) or `--slop`.
17. Restore `@U@` -> `_` in every output when a clean bank was used (`:585-612`).
18. Remove split parts unless `KeepSplits` (`:629-636`); remove the clean bank unless `KeepSplits` (`:639-646`); `cleanup_query` (`:649`).

## `sear2k` differences in the trace

- Round loop: `SKIP_MINUS` defaults 0, so minus rounds run until a round finds nothing (`sear2k:186-236`); `-m` skips them (`:221-223`).
- No empty-output guard: `bedtools merge` is fed directly (`sear2k:190-199`); no `set -e`.
- Final BED uses the bank `.fai` as genome file (`sear2k:261`), not a stripped one.
- Best mode dedupes by `(chrom,start,end,strand)` after sorting by bitscore (`sear2k:282-283`) rather than merging.

## Threading / parallelism model

Single process. Parts are processed in a serial `for` loop; the only
parallelism is `ssearch36 -T` (and `minimap2 -t`, `mafft --thread`), all set
from `THREADS` or `nproc`. No `xargs -P`, no background jobs. Memory is
negligible (streams through awk/bedtools).

## Error handling conventions

`sear`: `set -euo pipefail`; most bedtools calls have `2>/dev/null`; the
best-mode block disables `pipefail` around `head` (`:511,531`); MAFFT empty
output is fatal (`:566-569`); `trap cleanup_query EXIT`. `sear2k`: no `set -e`,
errors propagate as empty files.
