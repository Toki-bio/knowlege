---
title: "ssearch36-fasta36 — digest"
type: reference
topic: Program digests
generated: true
---

> Type: reference · Claude Code skill `ssearch36-fasta36` · [source on GitHub](https://github.com/Toki-bio/knowlege/tree/main/.claude/skills/ssearch36-fasta36)

**When to use:** Internals of W. R. Pearson's FASTA36 package (ssearch36 full Smith-Waterman, fasta36 ktup heuristic, lalign36 Waterman-Eggert, ggsearch36/glsearch36). Use when reading, modifying, debugging, or explaining ssearch36/fasta36 source, its -f/-g/-r/-s/-z/-E/-N/-m/-T options, E-value statistics, how long chromosomes are chunked, how many HSPs per library sequence are reported, SIMD performance, or when tuning ssearch36/fasta36 for sensitive SINE / repeat searches of DNA consensus sequences against whole genomes.

---

{% raw %}
# ssearch36 / fasta36 internals

> Type: reference (facts verified against the pinned commit) with marked `notes`.

## Pinned source

| Field | Value |
|---|---|
| Repo | https://github.com/wrpearson/fasta36 |
| Version / tag | `36.3.8j1 Aug, 2026` (`src/initfa.c:538`); readme calls it v36.4.8j1 (`doc/readme.v36:17`) |
| Commit | `6d6718a` (full SHA in references/architecture.md) |
| Digested on | 2026-10-11 |
| Language / build | C99 + SSE2 intrinsics (via SIMDe on ARM), hand-written Makefiles (`make/Makefile.linux64_sse2`) |
| License | Apache 2.0 (`LICENSE`); Farrar SSE2 code BSD (`COPYRIGHT.sse2`) |

## What it does in three sentences

One shared driver (`comp_lib9.c`) reads a query and streams a library through
worker threads; each program is the driver linked against a different
`drop*.c` "comparison function" (`dropgsw2.c` = Smith-Waterman for ssearch36 and
lalign36, `dropnfa.c` = ktup heuristic for fasta36, `dropnnw2.c` = global for
gg/glsearch36). The scan phase returns one best score per library entry (per
strand), E-values are fitted to the whole score distribution (regression of
score vs log length by default), and only for entries below `-E` are alignments
computed, recursively finding further non-overlapping local alignments inside
the same entry. Library sequences longer than 150,000 residues are silently split
into 150-residue-overlapping chunks that are scored and reported as separate
entries.

## Where to look (read on demand)

| Question | File |
|---|---|
| Directory map, which .c builds which program, compile defines | `references/architecture.md` |
| pstruct, mngmsg, beststr, a_res_str, seq_record, getlib_str | `references/data-structures.md` |
| Numbered trace of an ssearch36 run; threading model | `references/control-flow.md` |
| Farrar striped SW, 8/16-bit lanes, overflow fallback, DNA matrix and gap defaults, -f/-g/-r/-s | `references/smith-waterman-simd.md` |
| E-values: -z modes, regression, shuffles, -Z, -E/-z interplay, behaviour with few long sequences | `references/statistics.md` |
| Chunking of huge sequences (-N), one score per chunk, how many HSPs per entry, -E "e1 e2", -d/-b/-K | `references/library-chunking-and-hsps.md` |
| -T, buffers, memory, expected speed, what is not SIMD | `references/threading-performance.md` |
| fasta36 ktup / -c / -Xo / -A for DNA, lalign36 (SIM, -K) | `references/fasta-heuristic-lalign.md` |
| -m 8 / 8C / 8CB / 9c / 10 columns, per-HSP coordinates | `references/io-formats.md` |
| Gotchas, docs that disagree with code | `references/pitfalls.md` |
| Build, run, test | `references/build-run-test.md` |

Start with `library-chunking-and-hsps.md` when the question is "does ssearch36
report every SINE copy in a chromosome", with `statistics.md` for E-values, and
with `architecture.md` when the task is unfamiliar.

## Ten facts worth knowing before touching the code

1. DNA defaults (set in `resetp()`, `src/initfa.c:1679-1747`): matrix `+5/-4`
   (`initfa.c:638-639,1733-1737`), gap open `-12`, extend `-4`
   (`initfa.c:1680-1686`), so a gap of length k costs 12+4k; Farrar's routine
   receives open = -(gdelval+ggapval) = 16 (`src/dropgsw2.c:795-800`). `-E`
   default drops from 10 to 2.0 for DNA and the secondary-alignment threshold
   `e_cut_r` becomes 0.001 (`initfa.c:1584-1587`).
2. Only SSE2 (and Altivec) SIMD exists; there is no AVX2 path (grep of `src/`
   and `make/`). ARM builds route the same SSE2 code through SIMDe
   (`src/smith_waterman_sse2.c:27-28`, `make/Makefile.linux64_simde_arm`).
3. The SIMD kernel returns only the score (`dropgsw2.c:829`). Alignments
   (coordinates, CIGAR) are recomputed later with scalar Gotoh code
   (`src/wm_align.c:37` `sw_walign`) for the hits below `-E` only.
4. 8-bit lanes first; if a chunk's score saturates at 255 it is redone in 16-bit,
   and after 500 chunks with more than 1/3 overflows the search switches to
   16-bit permanently (`dropgsw2.c:801-826`, `smith_waterman_sse2.c:421-424`).
5. Library entries longer than `MAXLIB` = 150,000 (`src/defs.h:93`) are split
   into chunks with a fixed 150-residue overlap (`comp_lib9.c:893-894`,
   `compacc2.c:155-174`); `-N` can only make chunks smaller. Each chunk is an
   independent entry labelled `name:k` (`mshowbest.c:493-496`) but coordinates
   are reported in whole-sequence space via `l_offset` (`mshowbest.c:633-636`).
6. The scan keeps one best score per chunk per strand (`compacc2.c:2490-2611`).
   Multiple HSPs appear only in the alignment phase: `nsw_malign()` recursively
   aligns the library segments left and right of each alignment while the
   sub-score beats `repeat_thresh` (`wm_align.c:202-296`), which is the raw
   score equivalent of `e_cut_r` (`comp_lib9.c:1534-1550`). No count cap in
   ssearch36; `-K` applies to lalign36 only (`initfa.c:1201-1204`).
7. Default statistics `-z 1`: weighted regression of best score on ln(length)
   over up to 60,000 sampled entries (`scaleswn.c:429-537`, `defs.h:106`). With
   fewer than 500 scored entries the program automatically shuffles every
   library entry and refits (`comp_lib9.c:1079-1096`); with fewer than 10 it
   falls back to Altschul-Gish parameters (`scaleswn.c:289`).
8. For DNA the E-value database size is total residues / library entry length
   unless `-Z` is given (`scaleswn.c:2694-2699`); `-Z` for DNA:DNA only started
   working in July 2025 (`doc/readme.v36:34-38`).
9. Work unit for threads is a buffer of ~5000 entries (`comp_lib9.c:821-830`),
   so a genome of 20,000 150-kb chunks gives at most a few buffers per thread;
   `-T` defaults to all CPUs (`doinit.c:255`).
10. fasta36 on DNA: ktup 6 is the default and is a *positional* third argument,
    not a flag (`initfa.c:1469-1476,1696`); DNA alignments are banded (width 16)
    unless `-A` (`initfa.c:1711-1715`, `dropnfa.c:1559-1600`); initial words are
    not x-drop-extended for DNA (`dropnfa.c:607`).

## Open questions / unverified

- Throughput (GCUPS) is not stated anywhere in the source or docs; only
  "10-20 fold" SIMD speed-up (`doc/fasta36.1:67-70`) and "roughly 50% faster"
  for 8-bit vs 16-bit (`dropgsw2.c:819`). Measure on your hardware.
- Whether fasta36 with ktup 3-4 recovers ssearch36 sensitivity on 20-30%
  diverged SINEs is a guide claim (`doc/fasta_guide.tex:125-128`), not tested here.
- `-m 10` field names (`; al_start` etc.) come from `c_dispn.c:155,176`; the
  guide says `-m 10` is "no longer reliably tested" (`fasta_guide.tex:991`).

## Change log

- 2026-10-11: initial digest at `6d6718a`.
{% endraw %}

## Reference pages

- [architecture](architecture/) — Architecture
- [build-run-test](build-run-test/) — Build, run, test
- [control-flow](control-flow/) — Control flow
- [data-structures](data-structures/) — Core data structures
- [fasta-heuristic-lalign](fasta-heuristic-lalign/) — fasta36 heuristic vs ssearch36; lalign36
- [io-formats](io-formats/) — Input and output formats
- [library-chunking-and-hsps](library-chunking-and-hsps/) — Library handling: huge sequences, scores per entry, multiple HSPs
- [pitfalls](pitfalls/) — Pitfalls and gotchas
- [smith-waterman-simd](smith-waterman-simd/) — Smith-Waterman implementation (ssearch36)
- [statistics](statistics/) — Statistics: how E-values are estimated
- [threading-performance](threading-performance/) — Threading, performance, memory
