---
title: "bwa: Seeding (SMEMs) and chaining in `bwa mem`"
type: reference
generated: true
---

> Part of the [bwa](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/bwa/references/seeding-chaining.md)

{% raw %}
# Seeding (SMEMs) and chaining in `bwa mem`

## 1. The SMEM search — `bwt_smem1a` (`bwt.c:289-351`)

A super-maximal exact match (SMEM) is an exact match between query and index that cannot be extended in either direction and is not contained in another such match. `bwt_smem1a(bwt, len, q, x, min_intv, max_intv, mem, tmp)` returns all SMEMs **covering query position x**, and the position after the longest forward match (so the caller can jump ahead).

1. Start with the 1-base interval of `q[x]` (`bwt.c:301`, `bwt_set_intv`).
2. Forward phase (`304-321`): extend right one base at a time with `bwt_extend(.., is_back=0)` (`310`), using the complement because the forward extension is done on the reverse-complement interval `x[1]`. Each time the interval size `x[2]` shrinks, push the previous interval to `curr` (`311-312`). Stop when the size would fall below `min_intv` (`313`), at an N (`316-319`), at the end (`321`), or — only if `max_intv > 0` — as soon as the interval is already smaller than `max_intv` (`305-307`). `curr` is reversed so the longest (smallest-interval) match comes first (`322`).
3. Backward phase (`326-345`): for each candidate (longest first) try extending left (`330`). A candidate that can no longer be extended (or hits the query start / an N / `min_intv`) is emitted **only if no longer candidate was already extended at this step** (`curr->n == 0`, `332`) and it is not contained in the previously emitted SMEM (`333`). Candidates that still extend are deduplicated by interval size (`338-341`).
4. Result: SMEMs sorted by start (`346`), each as `bwtintv_t` with `x[2]` = number of occurrences on both strands and `info = start<<32 | end`.

Key consequences:
- Only *maximal* matches are reported. If the consensus shares a 19-mer with 10,000 copies but a 31-mer with 3 of them, the SMEM at that locus is the 31-mer with `x[2] = 3`. The 10,000-copy 19-mer is invisible to pass 1.
- `min_intv` sets a floor on occurrences (passing `x[2]+1` of a parent SMEM is how re-seeding asks for "shorter but more frequent").
- `bwt_smem1` is `bwt_smem1a` with `max_intv = 0` (`353-356`); the comment at `288` says `max_intv` is not used by BWA-MEM (the LAST-like pass uses `bwt_seed_strategy1` instead).

`bwt_seed_strategy1` (`358-379`): forward-only extension from `x`; returns the first match with length `> min_len` (`i - x >= min_len`) **and** `x[2] < max_intv` (`370`). It yields at most one seed per start and does not look for maximality.

## 2. The three passes — `mem_collect_intv` (`bwamem.c:140-188`)

| Pass | Lines | What it emits | Filter |
|---|---|---|---|
| 1 SMEM | `147-157` | all SMEMs with `min_intv = 1`, advancing `x` to the end of the longest match each time (`149`) | keep if length `>= -k` (`153`) |
| 2 re-seed | `159-168` | for each pass-1 SMEM with length `>= split_len = round(k * -r)` (`144`; 29 for k=19) **and** occurrences `<= -s split_width` (10): `bwt_smem1` from the SMEM midpoint with `min_intv = x[2]+1` (`164`) | length `>= k` (`166`). One level only — children are not re-seeded. Parents with > 10 occurrences are *never* re-seeded (`163`). |
| 3 LAST-like | `170-185` | if `-y max_mem_intv > 0`: tile the query with `bwt_seed_strategy1(.., min_len=k, max_intv=-y)` seeds (`176`), each `x` advancing to the end of the seed found | by construction every seed has `< -y` (20) occurrences and length `> k` |

Then sorted by `info` (`187`). Note: all three passes are biased toward *low-occurrence* seeds. There is no pass that systematically produces high-copy k-mers, except pass 3 with a very large `-y` (then `bwt_seed_strategy1` stops at the first length `> k` and emits a `(k+1)`-mer with whatever occurrence count it has).

## 3. From seeds to chains — `mem_chain` (`bwamem.c:277-341`)

1. Query shorter than `-k` -> no chains (`286`).
2. `frac_rep` (`291-298`): fraction of the query covered by seeds whose `x[2] > -c max_occ`; stored on every chain (`336`) and later multiplies MAPQ by `1 - frac_rep` (`1004`). For a consensus where every seed is high-copy, `frac_rep -> 1` and MAPQ -> 0.
3. **Occurrence cap (`-c`) is a subsample, not a drop** (`304-305`):
   ```
   step = x[2] > max_occ ? x[2] / max_occ : 1;
   for (k = count = 0; k < x[2] && count < max_occ; k += step, ++count)
   ```
   A seed with 10,000 copies contributes exactly 500 positions: every 20th entry of its SA interval. SA order is lexicographic by the suffix following the seed, so the chosen 500 are effectively arbitrary genomic copies; a different seed of the same family picks a different 500. The usage text "skip seeds with more than INT occurrences" (`fastmap.c:278`) and the man page "Discard a MEM" (`bwa.1:203-205`) describe older behaviour; the dropped line is still visible, commented out, at `303`.
4. Each sampled occurrence: `bwt_sa` locate (`309`, up to 31 LF steps), `bns_intv2rid` (`312`; seeds crossing a sequence boundary or the strand boundary are discarded, `313`), then look up the nearest chain by reference start in a B-tree (`315`) and `test_and_merge` (`216-237`): same `rid`, same strand (`225`), not fully contained in an existing seed (`223-224`), reference offset `y >= 0`, `|x - y| <= -w` (diagonal within the band), and both gaps `< -G max_chain_gap` (`228`). Otherwise a new chain (`318-325`).

Chain weight `mem_chain_weight` (`239-258`) = min over query/reference of the number of bases covered by seeds (overlaps counted once).

## 4. Chain filtering — `mem_chain_flt` (`bwamem.c:353-411`)

1. Drop chains with weight `< -W min_chain_weight` (default 0 -> none) (`359-365`).
2. Sort by weight descending (`367`). The heaviest is always kept (`369`).
3. For each other chain i, compare with every kept chain j (`371-393`). "Significant overlap" = query overlap `>= -X mask_level (0.5) * min(len_i, len_j)` and `min_l < max_chain_gap` (`381`); the ALT rule at `377`. On significant overlap: remember the first shadowed chain for MAPQ (`383`); **drop i only if `w_i < w_j * -D drop_ratio` AND `w_j - w_i >= 2k`** (`384-385`). Else i is kept with `kept = 2` (overlapping) or `3` (clean) (`389-392`).
4. The first chain shadowed by each kept chain is also retained (`kept = 1`, `394-397`) so `sub` can be estimated.
5. At most `-N max_chain_extend` chains with `kept < 3` survive (`399-404`); default unlimited.

So for a repeat query, every copy whose chain weighs at least half of the best chain's weight (or within 2k of it) is kept and extended. With `-D 0` nothing is dropped for weight. Chains are not limited by count unless `-N` is set.

`mem_flt_chained_seeds` (`624-641`) runs seed-level SW only when `5.5*ln(l) <= 0.05*l` (queries longer than ~1 kb); irrelevant for 100-500 bp queries.

## 5. Why a 300 bp query at 75% identity rarely has a 19-mer exact seed (math)

Model: per-base identity p between consensus and a given copy, independent sites (indels count as mismatches). Expected number of maximal exact runs of length >= k in n sites ~ `(n-k)(1-p)p^k + p^k`; P(at least one) ~ `1 - exp(-E)`.

| p (identity) | k | p^k | E (n=300) | P(>= 1 seed in that copy) |
|---|---|---|---|---|
| 0.90 | 19 | 0.135 | 3.8 | 98% |
| 0.85 | 19 | 0.046 | 1.9 | 85% |
| 0.80 | 19 | 0.0144 | 0.82 | 56% |
| 0.75 | 19 | 0.0042 | 0.30 | 26% |
| 0.70 | 19 | 0.0011 | 0.10 | 9% |
| 0.75 | 15 | 0.0134 | 0.95 | 61% |
| 0.75 | 13 | 0.0238 | 1.7 | 82% |
| 0.75 | 11 | 0.0422 | 3.0 | 95% |

And a seed is only the first hurdle: it must also (a) be reported by one of the three passes (pass 1 reports the *longest* match only; pass 2 only below 10 occurrences; pass 3 only below `-y` 20 occurrences), (b) survive the `-c` subsample, (c) form a chain that is >= 50% of the best chain, and (d) extend to score >= `-T`.

Random background for a k-mer in a 3.1 Gbp genome (N = 6.2e9 indexed positions, 4^k k-mers): expected occurrences `N/4^k`: k=19 -> 0.02; k=15 -> 5.8; k=13 -> 92; k=12 -> 370; **k=11 -> 1,480 (> 500)**; k=10 -> 5,900. With `-k 11` on a mammalian genome essentially every seed exceeds `-c 500`, so every seed is subsampled, `frac_rep -> 1`, MAPQ -> 0 for all hits, and each query spawns ~`(#seeds) x 500` chains, mostly single-seed noise chains of weight 11 that the filter keeps whenever the best chain is weak (`w_best * 0.5 <= 11` or `w_best - 11 < 22`). Each surviving chain costs one `ksw_extend2` left and right. Correctness is protected by `-T` (random 11-mer extensions seldom reach 30), but run time grows by orders of magnitude and output MAPQ carries no information. `-k` also changes the chain-filter slack `2k` (`384`), the default `sub` in MAPQ (`984`) and the re-seed length (`144`).

`notes:` `-k 11..13` is about the floor where exact seeds exist for 75%-identity copies; it is only usable together with `-c` raised to the family copy number (otherwise the subsample decides which copies you see) and, on gigabase genomes, with the expectation of seconds per query.

## 6. What happens to one seed hitting 10,000 genomic copies

1. If any copy shares a longer exact match, pass 1 reports that longer match (few copies) instead; the 10,000-copy k-mer only appears if it is itself maximal.
2. Pass 2 skips it (`x[2] > 10`, `163`); pass 3 skips it (`x[2] >= 20`, `bwt.c:370`).
3. If it is reported: `frac_rep` grows by its length (`294-297`); 500 of 10,000 positions are located (`304-309`); each becomes a 1-seed chain of weight k (unless another sampled seed lands on the same diagonal).
4. 9,500 copies get no seed from it. Unless another seed happens to sample them, they are never chained, never extended, never output — regardless of `-a`, `-T`, `-D`, `-h`.
{% endraw %}
