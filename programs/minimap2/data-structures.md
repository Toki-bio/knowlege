---
title: "minimap2: Core data structures"
type: reference
generated: true
---

> Part of the [minimap2](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/minimap2/references/data-structures.md)

{% raw %}
# Core data structures

All line numbers at commit `3c28777`.

## `mm128_t` — the universal 128-bit record (`minimap.h:77`)

`{uint64_t x, y}`. Its meaning depends on stage:

| Stage | `x` | `y` | Where set |
|---|---|---|---|
| Sketch output (query or ref) | `hash64(kmer) << 8 \| kmer_span` | `rid << 32 \| last_pos << 1 \| strand` | `sketch.c:112-113` |
| Index bucket entry (before `worker_post`) | same as sketch | same | `index.c:306-313` |
| Anchor `a[i]` (after seed lookup) | `rev << 63 \| rid << 32 \| ref_last_pos` | `flags << 40 \| q_span << 32 \| q_last_pos` (`seg_id` in bits 48-55) | `map.c:184-198`; layout comment `lchain.c:140-146` |

Anchor flag bits (`mmpriv.h:19-25`): `MM_SEED_LONG_JOIN` (1<<40),
`MM_SEED_IGNORE` (1<<41, set by `mm_filter_bad_seeds*`), `MM_SEED_TANDEM`
(1<<42, same minimizer as neighbour), `MM_SEED_SELF` (1<<43, all-vs-all),
`MM_SEED_SEG_SHIFT` 48. Reverse-strand anchors store the query position
flipped onto the reverse complement (`map.c:190`), so chaining treats both
strands as "forward".

Invariant: the anchor array is radix-sorted by `x` (`map.c:202`), i.e. by
strand, then rid, then ref position; chaining relies on this (`lchain.c:172`).

## `mm_idx_t` (`minimap.h:88-100`) and `mm_idx_bucket_t` (`index.c:28-33`)

| Field | Meaning | Invariant / note |
|---|---|---|
| `b, w, k, flag` | bucket bits (14 default, capped at 2k `index.c:53`), window, k-mer, `MM_I_HPC\|MM_I_NO_SEQ\|MM_I_NO_NAME` | Baked into `.mmi`; override CLI on load `index.c:626-627` |
| `n_seq, seq[]` | `mm_idx_seq_t {name, offset, len, is_alt}` | `offset` into `S` |
| `S` | 4-bit packed concatenated reference | absent when `MM_I_NO_SEQ` (`index.c:333,356-361`) |
| `B[1<<b]` | buckets; minimizer hash's low `b` bits pick the bucket (`index.c:97,310`) | |
| `B[i].a` | temp `mm128_v` of all minimizers in the bucket during build | freed after `worker_post` (`index.c:276-277`) |
| `B[i].p` | `uint64_t[]` of positions (`y` values) for minimizers occurring >1 time, grouped and sorted by position | `index.c:263-266` |
| `B[i].h` | khash `minimizer>>b<<1 -> value`; key low bit 1 = singleton with `y` stored inline in the value; else value = `start_p << 32 \| count` | `index.c:19-21, 256-266`; lookup `mm_idx_get` `index.c:93-110` |
| `I, spsc, J` | BED intervals / splice scores / junction jumps (splice modes only) | |
| `km, h` | arena; name->id hash (`mm_idx_index_name`) | |

Created by: `mm_idx_init` `index.c:50` via `mm_idx_gen` `index.c:389` or
`mm_idx_load` `index.c:516`. Freed by `mm_idx_destroy` `index.c:62`.
Memory per element: 16 B per minimizer while in `B[i].a`, then 8 B per
occurrence in `p[]` for non-singletons plus one khash slot (16 B key+val, load
factor 0.77 `khash.h:180`) per distinct minimizer. `mm_idx_stat`
(`index.c:112-134`) prints distinct count, singleton %, mean occurrence and
mean spacing (~`(w+1)/2`).

## `mm_seed_t` (`mmpriv.h:45-51`)

One query minimizer that was found in the index.

| Field | Meaning | Set at |
|---|---|---|
| `n` | number of reference occurrences | `seed.c:45` |
| `q_pos` | `last_pos << 1 \| strand` | `seed.c:40,45` |
| `q_span:31, flt:1` | k-mer span (differs from k with HPC); `flt=1` -> dropped as too frequent | `seed.c:45`, `seed.c:86-91,109-111` |
| `seg_id:31, is_tandem:1` | fragment segment; identical to neighbouring minimizer | `seed.c:45-48` |
| `cr` | pointer into `B[].p` or the inline singleton value | `seed.c:42-45` |

Lifetime: inside `mm_collect_matches` (`seed.c:98-132`); compacted to kept
seeds, then consumed by `collect_seed_hits[_heap]` (`map.c:102-204`) and freed.
Side outputs: `rep_len` (total query bp covered by filtered seeds) and
`mini_pos[]` (span<<32|pos of kept seeds, used by `mm_est_err`).

## Chaining outputs: `u[]` and compacted `a[]` (`lchain.c:140-146`)

`u[i] = score << 32 | n_anchors`; `a[]` is reordered so chain i's anchors are
contiguous, chains ordered by first anchor's `x` (`compact_a` `lchain.c:78-111`).
Score is the DP score of the chain end minus the score at the backtrack stop
(`lchain.c:51,66`), and a chain is emitted only when `score >= min_sc &&
n_anchors >= min_cnt` (`lchain.c:52,67`). Backtracking stops when the score
drops by more than `max_drop = bw` from its running peak (`mg_chain_bk_end`
`lchain.c:9-25`, `lchain.c:151`).

## `mm_reg1_t` — a hit (`minimap.h:112-127`) and `mm_extra_t` (`minimap.h:103-110`)

| Field | Meaning | Set by |
|---|---|---|
| `id, parent` | `parent == id` -> primary; else index of the primary it is secondary to | `mm_gen_regs` `hit.c:77-78`, `mm_set_parent` `hit.c:125-186` |
| `cnt, as` | number of anchors and offset into `a[]` | `hit.c:81-82` |
| `score, score0` | chaining score (s1 tag); `score0` = original before splits | `hit.c:79` |
| `subsc, n_sub` | best secondary's chaining score (s2 tag); number of near-equal secondaries | `hit.c:169-177` |
| `qs, qe, rs, re, rev, rid` | coordinates, from first/last anchor (chain) or from DP (aligned) | `mm_reg_set_coor` `hit.c:23-38`; `align.c:893-895` |
| `mlen, blen` | "fuzzy" matched bases and block length from anchors (cols 10/11 of PAF before alignment; recomputed after) | `mm_cal_fuzzy_len` `hit.c:8-21`; `mm_update_extra` `align.c:254` |
| `mapq` | 0-60 | `mm_set_mapq2` `hit.c:432-485` |
| `div` | minimizer-based divergence estimate (dv tag) | `mm_est_err` `esterr.c:62` |
| `split, split_inv, inv, strand_retained, is_alt, sam_pri, ...` | bit flags | `hit.c:122`, `align.c:866`, `hit.c:269`, `hit.c:226` |
| `hash` | per-hit hash used to break score ties deterministically | `hit.c:64-65,80` |
| `p` (`mm_extra_t*`) | only with `-c/-a`: `dp_score` (AS), `dp_max` (ms, peak segment score used by `-s`), `dp_max2` (best overlapping alt), `n_ambi`, `cigar[]` | `align.c:794,861,869,886`; `mm_update_extra` |

Created by `mm_gen_regs` (calloc), extended by `mm_insert_reg`
(`align.c:973`), freed in the output step (`map.c:629`).

## `mm_mapopt_t` (`minimap.h:136-192`) and `mm_idxopt_t` (`minimap.h:130-134`)

Index options: `k, w, flag, bucket_bits, mini_batch_size (50 Mbp), batch_size
(-I, 8 Gbp)`. Mapping options of interest (full defaults in
`parameters-and-presets.md`): `mid_occ_frac, min_mid_occ, max_mid_occ, mid_occ,
max_occ, max_max_occ, occ_dist, q_occ_frac` (seed filtering); `bw, bw_long,
max_gap, max_gap_ref, max_chain_skip, max_chain_iter, min_cnt,
min_chain_score, chain_gap_scale, chain_skip_scale, rmq_*` (chaining);
`mask_level, mask_len, pri_ratio, best_n` (primary/secondary); `a, b, q, e, q2,
e2, transition, sc_ambi, zdrop, zdrop_inv, end_bonus, min_dp_max, min_ksw_len,
max_clip_ratio, max_sw_mat` (alignment); `mini_batch_size (-K), cap_kalloc`
(resources). `mid_occ` is the single derived value most behaviour depends on
(`options.c:73-79`).

## `mm_tbuf_t` (`minimap.h:207-212`)

Per-thread arena (`km`) plus `rep_len` and `frag_gap` outputs of the last
query. One per mapping thread per batch (`map.c:557-559`), destroyed in the
output step (`map.c:582`). Arena is reset if the largest block exceeds 256 MB
or capacity exceeds `--cap-kalloc` (500 MB) (`map.c:371-376`).

## Pipeline records (`map.c:403-423`)

`pipeline_t` (shared: options, index, file handles, `n_processed`) and
`step_t` (one mini-batch: `seq[]`, `n_seg/seg_off` fragment grouping,
`reg[][]`, `n_reg[]`, `rep_len[]`, `frag_gap[]`, per-thread `buf[]`). A batch
holds `-K` bases of query (`map.c:551-552`; `bseq.c:99-100`) and all of its
hits until written.
{% endraw %}
