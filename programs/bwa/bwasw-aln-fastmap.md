---
title: "bwa: `bwasw`, `bwa aln`, and `fastmap` for diverged ~300 bp queries"
type: reference
generated: true
---

> Part of the [bwa](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/bwa/references/bwasw-aln-fastmap.md)

{% raw %}
# `bwasw`, `bwa aln`, and `fastmap` for diverged ~300 bp queries

## `bwa bwasw` (BWA-SW, deprecated: `main.c:70`)

Algorithm (`bsw2_aln1_core`, `bwtsw2_aux.c:248-295`; `bsw2_core`, `bwtsw2_core.c:449-619`):

1. The read is turned into its own small FM-index (`bwtl_seq2bwtl`, `bwt_lite.c:13`), i.e. a DAWG of the read's substrings (connectivity hash, `bwtsw2_core.c:99,459`).
2. DP runs between read-DAWG nodes and genome prefix-trie nodes (SA intervals of the genome FM-index). For each read node `v` and each child base `tj` (`499`), every cell's affine-gap score is computed (`fill_cell`, `421-433`), and the genome trie is descended with `bwt_2occ4` (`551-563`). This is "seed-free" in the sense that there is no fixed k-mer: matching starts from length-1 intervals and grows.
3. Heuristics that make it tractable, all of which cut sensitivity to multi-copy families:
   - **Z-best** (`-z`, default 1): per read node a heap keeps only the `z` best-scoring genome-trie cells (`467-468,514,544-547`); `cut_tail(w, opt->z, ..)` discards the rest (`134`, `583`, `598`). With `z = 1` only the single best branch survives at each step; the man page says "Higher -z increases accuracy at the cost of speed" (`bwa.1:604-605`), and suggests `-z 10` for noisy PacBio (`bwtsw2_main.c:72-73`).
   - **Band** `-w` (default 50, code `bwtsw2_aux.c:58`; man page says 33, stale): cells with `|tlen - qlen| > bw` are killed (`491-495`); `update_opt` further clamps it by query length (`552-556`).
   - **Score threshold** `-T` (default 30 x a, code `57`; man page says 37) and `-c` coefficient 5.5: the per-read threshold is `max(T, c * ln(l))` (`550`); for l = 300, `5.5 * 5.7 = 31`. Hits are saved when `G >= t` (`save_hits`, `223-246`).
   - **Narrow hits only** (`-s`, `is`, default 3): hits whose genome SA interval is `<= 3` are saved separately (`save_narrow_hits`, `248-271`, `254`) and are the only ones chained (`bsw2_chain_filter`, `281`, comment "only unique seeds are chained") and left-extended (`283`). Wider intervals are collapsed to one representative per interval, flagged as a "random hit" (`bsw2_resolve_duphits`, `273-310`, `300-307`) and printed with MAPQ 0 and CIGAR `*` (`bwtsw2_aux.c:503`). So an exact substring present in > 3 copies is explicitly treated as repetitive and not located.
4. Extension with scalar `ksw_extend` around the hits (`133,161`), merge, overlap resolution by `mask_level` (`294`), MAPQ up to 250 (`423-429`), reverse strand searched only if the forward strand has few supporting seeds (`599-602`).
5. Scoring defaults `-a1 -b3 -q5 -r2` (`54-57`), i.e. B/A = 3 -> implied divergence ~1.2% (same formula as mem). Expected score/base at 75% identity: `0.75 - 0.75 = 0`.

Verdict for diverged multi-copy queries: `bwasw` is **less** suitable than `bwa mem`. It is tuned to find the unique best hit of a long, nearly-identical read; copies beyond the Z-best branch are pruned inside the DP, repetitive intervals (> `-s 3` copies) are deliberately not resolved to coordinates, and it is single-hit-per-region by design. Raising `-z` and `-s` and lowering `-b` helps marginally at a steep speed cost; the code is unmaintained ("DEPRECATED").

## `bwa aln` / `samse` (BWA-backtrack)

- Bounded-error search over the FM-index (`bwt_match_gap`, `bwtgap.c:109`) with at most `max_diff` differences, where `max_diff` is derived from the read length and `-n 0.04` (missing probability under 2% error, `bwa_cal_maxdiff`, `bwtaln.c:42-54,111`; defaults `bwtaln.c:24-40`): for 100 bp reads this is ~5 mismatches. Gap opens `-o 1`, extensions `-e 6`, seed (`-l 32`) allows `-k 2` differences. The priority stack is capped at `-m 2,000,000` entries (`bwtgap.c:139`) and search stops at `-R 30` equally best hits (`max_top2`). The man page itself says "I would not recommend to use BWA on data with r > 0.02" (`bwa.1:780-781`).
- For a 300 bp query at 25% divergence you would need ~75+ differences; the algorithm is exponential in `max_diff` and `-n` is capped by what the stack allows. Not usable.
- `samse -n N` reports extra hits in XA only if the total occurrences are `<= N+1`, else none at all (`bwase.c:57-58`).

## `bwa fastmap` as a cheap exact pre-filter

- Prints every SMEM (pass 1 only, no re-seeding) with its occurrence count, and genomic coordinates when the count is `<= -w` (default 20) (`fastmap.c:456-473`). `-l` minimum SMEM length (17), `-i` minimum interval size, `-I` maximum interval size (stop extending once a match is longer than `-l` and rarer than `-I`, `fastmap.c:437`, via `bwt_smem1a` `max_intv`), `-L` maximum MEM length.
- Memory: BWT + SA + annotation only (`444`), ~4.65 GB for human; no SW, so it is seeding cost only.
- What it gives you: for each query, the longest exact matches and **how many copies share each**. With `-w` raised to the family size (e.g. `-w 100000`) and `-l 11..15`, it will list the coordinates of every copy that shares a maximal exact match — but still only *maximal* matches (same SMEM caveat as mem: a 19-mer shared by 10,000 copies is hidden behind any longer match to fewer copies), and locating each occurrence costs a `bwt_sa` walk (`466`).
- Useful as: a survey of which consensus segments are exactly conserved and how repetitive they are (`EM` lines with large counts), a sanity check of what `bwa mem` could possibly seed, and a cheap way to pick `-k`/`-c`. Not a substitute for an alignment: no mismatches, no scores.
{% endraw %}
