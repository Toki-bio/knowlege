---
title: "bwa: Extension, scoring, dedup, primary/secondary/supplementary, MAPQ, output limits"
type: reference
generated: true
---

> Part of the [bwa](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/bwa/references/extension-scoring-output.md)

{% raw %}
# Extension, scoring, dedup, primary/secondary/supplementary, MAPQ, output limits

## 1. Extension kernel — `ksw_extend2` (`ksw.c:416-515`)

Scalar (no SIMD) banded affine-gap SW that extends from a seed with initial score `h0` (seed length x A, `bwamem.c:757`).

- Band `w` is first clamped to the longest gap the query could afford: `(qlen*max_match + end_bonus - o)/e + 1` (`434-443`).
- Per target row: band `[i-w, i+w+1)` (`452-453`), H/E/F recurrences with M separated so a cell cannot be reached through gap-into-gap (`469`); the row's max and its column are tracked (`473-474`).
- `gscore`: best score among cells in the **last query column** (`486-489`) = "reached the end of the query".
- Stop if the row max is 0 (`490`). **Z-drop** (`494-500`): if the row max is more than `zdrop` below the global max after discounting the diagonal drift `(|di - dj|) * e`, stop. With `-d 100 -A 1`, a run of net ~100 score loss ends the extension.
- The live band shrinks to the non-zero region each row (`502-505`), so cost is far below `qlen * tlen` for divergent sequence.
- Returns `max`, `qle/tle` (position of max), `gtle/gscore` (to-end), `max_off` (max diagonal offset, `491-493`).

In `mem_chain2aln` (`bwamem.c:658-812`): the reference window is the chain's seeds +/- the affordable gap (`668-683`, `cal_max_gap` `647-654`, capped at `2w`). Seeds are extended longest-first (`688-695`); a seed fully inside an existing region on the same diagonal is skipped (`697-732`). Left and right extension each try band `w` then `2w` if `max_off >= 0.75 w` (`749-760`, `777-788`, `MAX_BAND_TRY 2`). The to-end alignment is preferred when `gscore > local_score - pen_clip` (`762-768`, `790-796`): clipping penalty `-L` is this margin, not a score deduction (`bwa.1:247-253`). `seedcov` = seed bases inside the region (`801-805`).

The SSE2 striped kernels `ksw_u8`/`ksw_i16` (`ksw.c:122,255`, via `ksw_align2` `379`) are used only by `mem_seed_sw` (`bwamem.c:619`, long queries) and mate rescue (`bwamem_pair.c:137`). The final CIGAR is produced by banded global `ksw_global2` (`ksw.c:540-642`) through `bwa_gen_cigar2` (`bwa.c:148-234`), re-run with a doubled band up to 3 times if the global score falls short of the local one (`bwamem.c:1144-1152`).

## 2. Scoring (`-A -B -O -E -L`, `-T`)

- Matrix: `A` on the diagonal, `-B` off, `-1` for N (`bwa_fill_scmat`, `bwa.c:136-145`). Gap of length k costs `O + k*E` (`bwa.1:240-245`).
- Theory (`bwamem.c:49-70`): a mismatch penalty `B` relative to match `A` encodes an assumed substitution rate `e = 0.75 * exp(-ln(4) * B/A)` (`bwa.1:234`).

| B/A | implied divergence e | expected score per aligned base at identity p = 0.90 / 0.80 / 0.75 / 0.70 |
|---|---|---|
| 4 (default) | 0.3% | +0.5 / 0.0 / -0.25 / -0.5 |
| 2 | 4.7% | +0.7 / +0.4 / +0.25 / +0.1 |
| 1 (`-x pacbio/ont2d`) | 18.75% | +0.8 / +0.6 / +0.5 / +0.4 |

Expected score/base = `p*A - (1-p)*B` (gaps extra). At default scoring a 75%-identity copy has negative drift: extension beyond the seed loses score, Z-drop fires, the region stays a few tens of bp, and `-T 30` is rarely reached. 300 bp at 75% identity with `-A1 -B1 -O1 -E1` gives ~+150 before gaps.

- `-T` is applied only at output (`bwamem.c:1048`; "only affecting output", `bwamem.h:63`). `-A` rescales unspecified `-T -d -B -O -E -L -U` (`fastmap.c:125-139`).
- Presets (`fastmap.c:330-359`): `-x pacbio` = `-k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0`; `-x ont2d` = `-k14 -W20 -r10 ...`; `-x intractg` = `-B9 -O16 -L5`.

## 3. Dedup and merge — `mem_sort_dedup_patch` (`bwamem.c:463-515`)

Sorted by reference end; two regions on the same `rid` within `max_chain_gap` are redundant if they overlap by `> 0.95` (`mask_level_redun`) of the shorter on **both** reference and query (`481`) — the lower-scoring one is deleted. Otherwise `mem_patch_reg` (`432-461`) may merge co-linear pieces into one alignment via a global DP when the merged score is >= 90% of the predicted sum. Then sort by score (`504`) and remove exact duplicates (`505-513`). Distinct genomic copies are never "redundant" (different `rb`), so all survive this step.

## 4. Primary vs secondary vs supplementary

`mem_mark_primary_se` (`547-584`) sorts by score (ties by ALT status and a hash of read id, `423-424,556`) and runs `mem_mark_primary_se_core` (`519-545`): walking down the list, an alignment that overlaps an already-kept alignment on the **query** by `>= mask_level (0.5) * min(len)` becomes `secondary = index of that kept alignment` (`530-543`); the kept one records `sub = first competitor score` (`535`) and increments `sub_n` for competitors within `max(A+B, O+E)` of its score (`522-524,536-537`). Non-overlapping alignments stay primary (`542`).

In `mem_reg2sam` (`1033-1079`), for each region in score order:
- skip if `score < T` (`1048`);
- skip secondaries unless `-a` (`1049`); skip a secondary if `score < drop_ratio * parent score` (`1050`) — so `-a` alone shows only secondaries scoring >= 50% of their primary; use `-D 0` to lift that;
- the first printed region is the primary; later **primary** regions (other query parts) get flag 0x800 supplementary, or 0x100 with `-M` (`1057-1058`), and their MAPQ is capped to the first's unless `-q/-5` (`1059-1060`);
- secondaries: flag 0x100 and MAPQ 0 (`mem_reg2aln`, `1136-1137`), `XS` suppressed (`1056`), SEQ/QUAL `*` (`897-898`), still carry `AS`, `NM`, `MD`, CIGAR with soft clips.

Consequence for a multi-copy repeat query: every copy aligning to the same query interval is secondary to the single best copy. Without `-a` the user sees one SAM line (MAPQ ~0) and possibly an XA tag; with `-a` they see every copy with score >= max(T, 0.5 x best) as 0x100 lines.

## 5. XA tag — `mem_gen_alt` (`bwamem_extra.c:124-172`)

Candidates: regions whose `secondary_all` parent exists and whose score `>= XA_drop_ratio (-z 0.8) * parent score` (`get_pri_idx`, `116-121`). Per primary, if the candidate count exceeds `-h max_XA_hits` (5; 200 when an ALT hit is present) **no XA is written at all** (`144`). `-a` skips XA generation entirely (`bwamem.c:1041-1042`). XA entry format `chr,±pos,CIGAR,NM;` (`146-161`); `-u` writes `XB` with `,score,mapq` appended.

So for a repeat with hundreds of near-equal copies, default output is one line, MAPQ 0, **no XA**, unless `-h` is raised to at least the copy number (and `-z` lowered to include more diverged copies).

## 6. Mapping quality — `mem_approx_mapq_se` (`bwamem.c:982-1006`)

```
sub   = max(a->sub or k*A if none, a->csub)
if sub >= score: return 0
identity = 1 - (l*A - score) / (A+B) / l        (l = longer of query/ref span)
tmp = (l < 50) ? 1 : ln(50)/ln(l);  tmp *= identity^2
mapq = 6.02 * (score - sub)/A * tmp^2
mapq -= 4.343 * ln(sub_n + 1)
clamp to [0,60]
mapq *= (1 - frac_rep)
```
MAPQ is a phred-scaled score *difference* to the best overlapping competitor, shrunk by identity, length, number of near-ties, and the repetitive-seed fraction. For a repeat consensus, `sub ~ score` (other copies) -> 0, and `frac_rep ~ 1` -> 0 anyway. MAPQ therefore says nothing about alignment quality for repeat queries; use `AS`/`NM`.

## 7. Supplementary split alignments

Different query intervals (overlap < 50%) each produce a primary-class alignment; e.g. a composite element whose halves match different loci yields a primary + a supplementary line linked by `SA:Z` (`939-958`), hard-clipped unless `-Y` (`844-845`). `-5` makes the leftmost query piece primary (`mem_reorder_primary5`, `1008-1030`).
{% endraw %}
