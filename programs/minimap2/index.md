---
title: "minimap2 — digest"
type: reference
topic: Program digests
generated: true
---

> Type: reference · Claude Code skill `minimap2` · [source on GitHub](https://github.com/Toki-bio/knowlege/tree/main/.claude/skills/minimap2)

**When to use:** Internals of minimap2 (Heng Li's minimizer-seed / chain / ksw2-extend aligner for long reads, short reads, cDNA and genome-to-genome). Use when reading, modifying, debugging, or explaining minimap2 source, presets, parameters (-k/-w/-f/-n/-m/-s/-N/-p/-P), PAF/SAM output, or when deciding whether and how to tune minimap2 for sensitive repeat / SINE / transposon searches against whole genomes.

---

{% raw %}
# minimap2 internals

> Type: reference (facts verified against the pinned commit) with marked `notes`.

## Pinned source

| Field | Value |
|---|---|
| Repo | https://github.com/lh3/minimap2 |
| Version / tag | 2.31-r1302 (`MM_VERSION`, `minimap.h:8`) |
| Commit | `3c28777` (full SHA in references/architecture.md) |
| Digested on | 2026-10-11 |
| Language / build | C (C99, SSE2/SSE4.1/NEON intrinsics), Makefile |
| License | MIT (`LICENSE.txt`) |

## What it does in three sentences

minimap2 sketches every reference sequence into (w,k)-minimizers stored in a
hash table (`index.c`), then for each query collects the reference positions
("anchors") of the query's minimizers, deliberately *dropping* minimizers that
are too frequent in the reference (`seed.c`, `-f`). Anchors are chained with a
sparse dynamic program (`lchain.c`) under gap/bandwidth/min-count/min-score
limits, chains are ranked into primary/secondary hits (`hit.c`) and, only when
CIGAR output is requested (`-c`/`-a`), base-level alignment is filled in by the
ksw2 SIMD extension aligner (`align.c`, `ksw2_*`). It is a *mapper*: every
default is tuned to find the one true origin of a read and to suppress
repeats, which is the opposite of what a repeat/SINE search wants.

## Where to look (read on demand)

| Question | File |
|---|---|
| How is the code organized, what calls what | `references/architecture.md` |
| The structs everything revolves around (mm_idx_t, mm128_t anchors, mm_seed_t, mm_reg1_t, mm_mapopt_t) | `references/data-structures.md` |
| Numbered trace: index build, one query through seed/chain/align/output, threading | `references/control-flow.md` |
| Every preset's exact values, every flag the SINE question touches, memory formula | `references/parameters-and-presets.md` |
| PAF/SAM columns and tags, .mmi index layout | `references/io-formats.md` |
| Sensitive repeat search: which minimizers survive, the k/w/identity math, flags, what breaks, honest verdict | `references/repeat-search-tuning.md` |
| Gotchas, surprising defaults, stale docs | `references/pitfalls.md` |
| Build, run, test, debug flags (`--print-seeds`, `--dbg-seed-occ`, `-v3`) | `references/build-run-test.md` |

Read only the file that answers the current question. Start with
`architecture.md` when the task is unfamiliar; start with
`repeat-search-tuning.md` when the question is "why does minimap2 miss my
SINEs".

## Ten facts worth knowing before touching the code

1. Defaults are `-k15 -w10` (`options.c:8`), `-f 2e-4`, `-n 3`, `-m 40`,
   `-r 500,20000`, `-g 5000`, `-p 0.8`, `-N 5`, `-A2 -B4 -O4,24 -E2,1`,
   `-z 400,200` (`options.c:14-67`). `-s` (min DP score) defaults to
   `min_chain_score * a` = **80**, not 40 (`options.c:52`); `sr` sets 40, `asm*`
   and `map-hifi` set 200.
2. `-f FLOAT` becomes an absolute occurrence cutoff `mid_occ` computed from the
   index (`mm_idx_cal_max_occ`, `index.c:198-220`), then clamped into
   `[-U lo, -U hi]` = `[10, 1000000]` (`options.c:73-79`). `-f 0` therefore
   means `mid_occ = 1,000,000`, not "no filter". `-f INT` (>= 1) sets `mid_occ`
   directly and skips the clamp (`main.c:325-331`).
3. A query minimizer with more than `mid_occ` reference hits is dropped
   (`seed.c:109-111`), except that per streak of high-occurrence seeds the code
   rescues the `round(streak_len / occ_dist)` *lowest*-occurrence ones (at most
   128, each <= `max_max_occ` 4095) (`seed.c:56-96`; `-e` sets `occ_dist`=500).
   For a 300 bp query made entirely of repeat seeds that is **one** rescued
   seed. Dropped seed span is reported as the `rl:i` PAF tag (`seed.c:117-122`).
4. Chaining score between anchors = `min(q_span, gap) - (0.8*0.01*k*|dd| +
   0.5*log2(|dd|+1))` (`lchain.c:113-138`, `map.c:273`); a chain survives only
   if score >= `-m` **and** anchor count >= `-n` (`lchain.c:52,67`).
5. Primary/secondary is decided before alignment: a chain overlapping a better
   chain on the query by > `-M 0.5` of the shorter one is secondary
   (`hit.c:165`); secondaries are kept only if score >= `-p` * primary (or
   within `2k`) and fewer than `-N` so far (`hit.c:265-267`). `-p 0` disables
   this filter entirely (`hit.c:257`); `-P` skips primary assignment altogether
   (`map.c:208,219`).
6. Base alignment happens only with `-c`/`-a` (`map.c:217`); it runs ksw2
   `extz2` (single affine) or `extd2` (dual affine) (`align.c:356-359`), left
   extension -> gap filling between seeds (only when both gaps >= `min_ksw_len`
   200, `align.c:810`) -> right extension, with Z-drop and inversion re-test
   (`align.c:843-868`). After alignment `mm_filter_regs` drops hits with
   `cnt < -n`, `mlen < -m`, `dp_max < -s` (`hit.c:301-320`).
7. Index memory: 16 B per minimizer while collecting (`index.c:306-313`), then
   8 B per occurrence of non-singleton minimizers + ~16-32 B per distinct
   minimizer in khash (`index.c:239-272`), + 0.5 B/base of 4-bit sequence
   unless `MM_I_NO_SEQ` (`index.c:356-361`). Minimizer density ~ 2/(w+1);
   `-w 1` keeps every k-mer. Sequence is *not* stored when neither `-d` nor
   `-c`/`-a` is given (`main.c:352-353`).
8. Threading: `kt_pipeline` with 3 steps (read / map / write) on 2 threads (3
   with `-2`) (`map.c:677-678`); step 1 fans out over fragments with
   work-stealing `kt_for` on `-t` threads (`map.c:576`, `kthread.c:30-72`).
   `-K` (500 Mbp of query per mini-batch, `options.c:56`) bounds memory: all
   hits of a batch are held until output (`map.c:578-641`).
9. `-x` presets are applied in a first pass before any other option
   (`main.c:150-163`), so later flags override preset values. `-k/-w/-H` are
   baked into a `.mmi` and silently override command-line values
   (`index.c:626-627`).
10. `--sr`/`-x sr` is a different code path: query-length chaining gap,
    rechain with a higher `max_occ` when no chain forms (`map.c:293-316`),
    longest-exact-stretch seeding plus ungapped-first alignment with
    `end_bonus` (`align.c:664-669, 696-703, 823-834`); it asserts no HPC
    (`align.c:655`). `--rmq` swaps the DP chainer for the minigraph RMQ
    chainer (`map.c:275-277`) and refuses to combine with `--sr`/`--splice`
    (`options.c:209-213`).

## Open questions / unverified

- Actual `mid_occ` values for a given genome at k=11..15 were not measured
  here; run with `-v3` to see `mid_occ = N` (`options.c:81-82`).
- Conserved-minimizer expectations in `repeat-search-tuning.md` use the
  2/(w+1) density approximation and an independent-substitution model; they
  were not validated by simulation in this digest.
- Wall-clock and RSS numbers for `-k11 -w1` on a 3 Gbp genome are derived from
  the per-item byte costs above, not measured.
- SAM/PAF writer details beyond the tag list (`format.c:397-458`) were skimmed.

## Change log

- 2026-10-11: initial digest at `3c28777` (2.31-r1302).
{% endraw %}

## Reference pages

- [architecture](architecture/) — Architecture
- [build-run-test](build-run-test/) — Build, run, test
- [control-flow](control-flow/) — Control flow
- [data-structures](data-structures/) — Core data structures
- [io-formats](io-formats/) — Input and output formats
- [parameters-and-presets](parameters-and-presets/) — Parameters and presets
- [pitfalls](pitfalls/) — Pitfalls and gotchas
- [repeat-search-tuning](repeat-search-tuning/) — Sensitive repeat / SINE search with minimap2
