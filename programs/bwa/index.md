---
title: "bwa — digest"
type: reference
topic: Program digests
generated: true
---

> Type: reference · Claude Code skill `bwa` · [source on GitHub](https://github.com/Toki-bio/knowlege/tree/main/.claude/skills/bwa)

**When to use:** Internals of BWA (lh3/bwa 0.7.19 — index, mem, fastmap, bwasw, aln/samse/sampe). FM-index layout and on-disk sizes, SMEM seeding, chaining, ksw extension, scoring, MAPQ/secondary/supplementary logic, threading. Use when reading, modifying, debugging or explaining bwa source, parameters, index files or SAM output, or when judging whether and how `bwa mem` can serve as a sensitive aligner for diverged or high-copy queries (e.g. a SINE/repeat consensus against a whole genome).

---

{% raw %}
# BWA internals

> Type: reference (facts verified against the pinned commit) with marked `notes`.

## Pinned source

| Field | Value |
|---|---|
| Repo | https://github.com/lh3/bwa |
| Version / tag | 0.7.19-r1273 (`main.c:33`); last release 22 Mar 2025 (`NEWS.md`) |
| Commit | `d82444c` (full SHA in references/architecture.md), "recommend bwa-mem3", 2026-08-07 |
| Digested on | 2026-10-11 |
| Language / build | C (C99 + pthreads + zlib), single `Makefile`, `-O3`; no test suite, CI only builds |
| License | GPLv3 for the package (`COPYING`); most source files carry MIT headers |

## What it does in three sentences

BWA builds a 2-bit packed reference plus an FM-index (BWT with 64-bit Occ
checkpoints every 128 bases and a suffix array sampled every 32) of the genome
and its reverse complement, costing about 1.75 bytes per genome base on disk
and in RAM (`references/index-and-memory.md`). `bwa mem` seeds a query with
super-maximal exact matches (SMEMs, `bwt.c:289`), subsamples each seed to at
most `-c 500` genomic occurrences, chains co-linear seeds, keeps only chains
whose weight is competitive with the best overlapping chain, extends each with
a scalar banded Smith-Waterman using Z-drop (`ksw.c:416`), and reports the best
alignment per query region as primary with everything overlapping it as
secondary (MAPQ 0, hidden unless `-a`). `bwasw` (deprecated) does DP between
the read's own BWT and the genome FM-index, and `bwa aln` does bounded-error
backtracking; both are tuned for near-identical reads.

## Where to look (read on demand)

| Question | File |
|---|---|
| Code layout, who calls whom, entry points | `references/architecture.md` |
| The structs everything revolves around (`bwt_t`, `mem_opt_t`, `mem_alnreg_t`...) | `references/data-structures.md` |
| Numbered trace of `index`, `mem`, `fastmap`, `bwasw`, `aln`; threading | `references/control-flow.md` |
| Index files `.bwt .sa .pac .ann .amb`, sizes and formulas, build algorithms, RAM | `references/index-and-memory.md` |
| SMEM algorithm, re-seeding, `-k -r -s -y -c`, chaining and chain filtering, seed-hit probability math | `references/seeding-chaining.md` |
| `ksw_extend2`, Z-drop, scoring `-A -B -O -E -L -T`, dedup, primary/secondary/supplementary, MAPQ, XA/`-h`/`-a` | `references/extension-scoring-output.md` |
| `bwasw`, `bwa aln/samse`, `fastmap` as a pre-filter | `references/bwasw-aln-fastmap.md` |
| Is `bwa mem` usable as a repeat finder; closest flag set; cost | `references/repeat-finding-assessment.md` (notes) |
| Input/output formats, SAM tags | `references/io-formats.md` |
| Gotchas, stale docs, surprising defaults | `references/pitfalls.md` |
| Build, run, debug | `references/build-run-test.md` |

Start with `architecture.md` for unfamiliar tasks; for the repeat/SINE
question read `seeding-chaining.md` then `repeat-finding-assessment.md`.

## Ten facts worth knowing before touching the code

1. Index size is linear: `.bwt` ~1.0 B/base, `.sa` ~0.5 B/base, `.pac` ~0.25 B/base
   (both strands are indexed, so `seq_len = 2*l_pac`). A 3.1 Gbp genome ~5.4 GB on
   disk and in RAM for `mem` (`bwtindex.c:74,155-156`, `bwt.c:73`, `bntseq.c:316`).
2. `bwa index` hard-codes SA sampling 32 (`bwtindex.c:316`); only the low-level
   `bwa bwt2sa -i` changes it (`bwtindex.c:191`). Occ interval 128 is a compile-time
   constant that must not change (`bwt.h:36-38`).
3. Seeding reports only *maximal* exact matches: a 19-mer shared by 10,000 copies is
   not emitted if any copy extends it further; re-seeding fires once, only for SMEMs
   with <= `-s 10` occurrences and length >= 1.5k (`bwamem.c:163-164`).
4. `-c 500` does **not** discard a seed (despite the usage text): it keeps at most 500
   occurrences, every `n/500`-th SA entry (`bwamem.c:304-305`), and such seeds raise
   `frac_rep`, which scales MAPQ down (`bwamem.c:291-298,1004`).
5. Chain filter drops a chain only if it overlaps a heavier chain by >= 50% of the
   shorter one *and* weighs < `-D 0.5` of it *and* is >= 2k lighter (`bwamem.c:381-385`).
6. Extension is **scalar** banded SW with Z-drop (`ksw_extend2`, `ksw.c:416-515`); the
   SSE2 striped `ksw_u8/ksw_i16` is used only by seed-SW and mate rescue.
7. Default scoring `-A1 -B4` models ~0.3% divergence (`e = .75*exp(-ln4*B/A)`,
   `bwamem.c:57`); expected score/base is 0 at 80% identity and negative at 75%.
8. All alignments overlapping the primary by >= 50% on the query become secondary:
   flag 0x100, MAPQ 0, SEQ `*`, printed only with `-a` and only if score >= `-T` and
   >= `-D`*primary (`bwamem.c:528-543,1048-1050,1136-1137`, `897-898`).
9. XA is all-or-nothing: if more than `-h 5` hits score >= 80% of the primary, no XA
   tag is written at all (`bwamem_extra.c:144`). `-a` disables XA (`bwamem.c:1041`).
10. Threads: a 2-stage `kt_pipeline` (read / align / write) wraps `kt_for` with work
    stealing over reads; batch = 10 Mbp x threads unless `-K` (`fastmap.c:394-395`,
    `kthread.c:25-61`).

## Open questions / unverified

- Absolute throughput (GCUPS, reads/s) is not derivable from the code; figures in
  `repeat-finding-assessment.md` are order-of-magnitude notes.
- `-a rb2` (rope-based BWT, `bwtindex.c:95-119`) memory and time not measured.
- `bwasw` RAM (whether `.pac` is memory-mapped or read) not traced beyond `bwtsw2_main.c:83`.
- Behaviour of bwa-mem2 / bwa-mem3 (recommended by README) not covered; same algorithm
  claimed, different index format.

## Change log

- 2026-10-11: initial digest at `d82444c`.
{% endraw %}

## Reference pages

- [architecture](architecture/) — Architecture
- [build-run-test](build-run-test/) — Build, run, test
- [bwasw-aln-fastmap](bwasw-aln-fastmap/) — `bwasw`, `bwa aln`, and `fastmap` for diverged ~300 bp queries
- [control-flow](control-flow/) — Control flow
- [data-structures](data-structures/) — Core data structures
- [extension-scoring-output](extension-scoring-output/) — Extension, scoring, dedup, primary/secondary/supplementary, MAPQ, output limits
- [index-and-memory](index-and-memory/) — Index construction, on-disk size and memory
- [io-formats](io-formats/) — Input and output formats
- [pitfalls](pitfalls/) — Pitfalls and gotchas
- [repeat-finding-assessment](repeat-finding-assessment/) — Is `bwa mem` usable as a sensitive repeat finder? (Type: notes)
- [seeding-chaining](seeding-chaining/) — Seeding (SMEMs) and chaining in `bwa mem`
