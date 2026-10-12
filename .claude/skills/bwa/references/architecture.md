# Architecture

Full commit SHA: `d82444c17edc2384420409f85557c6ae84019732` (https://github.com/lh3/bwa, 2026-08-07, "recommend bwa-mem3"). Version string `0.7.19-r1273` (`main.c:33`).

## Directory map (flat repository, ~17k lines of C)

| Group | Files | Role |
|---|---|---|
| Driver | `main.c` | Sub-command dispatch (`main.c:98-114`); prints `@PG` line, timing |
| Index construction | `bwtindex.c`, `bwt_gen.c` (BWT-SW incremental construction, from Wong/Lam), `is.c` (SA-IS), `rope.c`/`rle.c` (rb2 rope BWT), `QSufSort.c` | `bwa index`, `fa2pac`, `pac2bwt`, `bwtupdate`, `bwt2sa` |
| Reference packing | `bntseq.c/.h` | FASTA -> `.pac/.ann/.amb`; coordinate mapping (`bns_depos`, `bns_pos2rid`, `bns_fetch_seq`) |
| FM-index | `bwt.c/.h` | Occ/LF, SA lookup, bidirectional extension, SMEM search (`bwt_smem1a`), third-round seeding (`bwt_seed_strategy1`), load/dump |
| Shared infra | `bwa.c/.h` (index loading, batch FASTQ reader, scoring matrix, CIGAR via global DP, SAM header), `bwashm.c` (POSIX shm index), `utils.c`, `kstring.c`, `kseq.h`, `kvec.h`, `ksort.h`, `kbtree.h`, `khash.h`, `malloc_wrap.c` | |
| Alignment kernels | `ksw.c/.h` | SSE2 striped local SW (`ksw_u8`, `ksw_i16`, `ksw_align2`), scalar banded extension with Z-drop (`ksw_extend2`), scalar banded global with CIGAR (`ksw_global2`). `neon_sse.h`/`scalar_sse.h` provide portability shims (`ksw.c:29-35`) |
| Threads | `kthread.c` | `kt_for` (work-stealing parallel for) and `kt_pipeline` (multi-step pipeline) |
| BWA-MEM | `bwamem.c` (core), `bwamem_extra.c` (SMEM iterator API, `mem_align1`, XA generation), `bwamem_pair.c` (insert-size estimation, mate rescue, pairing), `fastmap.c` (`bwa mem` and `bwa fastmap` CLIs) | |
| BWA-SW | `bwtsw2_main.c` (CLI), `bwtsw2_core.c` (DP over read-BWT x genome FM-index), `bwtsw2_aux.c` (per-read driver, extension, SAM), `bwtsw2_chain.c`, `bwtsw2_pair.c`, `bwt_lite.c` (BWT of the read) | |
| BWA-backtrack | `bwtaln.c` (`aln`), `bwtgap.c` (bounded-error backtracking search), `bwase.c` (`samse`), `bwape.c` (`sampe`), `bwaseqio.c`, `bamlite.c` | |
| Misc | `pemerge.c`, `maxk.c`, `example.c` (libbwa API example), `bwakit/` (scripts), `xa2multi.pl`, `qualfa2fq.pl` | |
| Docs | `README.md`, `README-alt.md`, `NEWS.md`, `ChangeLog`, `bwa.1` (man page) | |

Build products: `libbwa.a` from `LOBJS` (`Makefile:7-8`) plus the `bwa` binary from `AOBJS` (`Makefile:9-13`). `make bwamem-lite` builds `example.c` against the library.

## Module dependency sketch (top-down)

```
main.c
 ├─ bwtindex.c ── bntseq.c (pack) ── bwt_gen.c | is.c | rope.c (BWT) ── bwt.c (occ update, SA sampling)
 ├─ fastmap.c:main_mem ── kthread.c:kt_pipeline ── bwa.c:bseq_read
 │       └─ bwamem.c:mem_process_seqs ── kthread.c:kt_for
 │             ├─ worker1: mem_align1_core ── mem_chain (bwt.c SMEM + bwt_sa) ── mem_chain_flt
 │             │        ── mem_chain2aln (ksw.c:ksw_extend2) ── mem_sort_dedup_patch (bwa.c:bwa_gen_cigar2)
 │             └─ worker2: mem_mark_primary_se ── mem_reg2sam ── mem_reg2aln (bwa_gen_cigar2 -> ksw_global2)
 │                                               └─ bwamem_extra.c:mem_gen_alt (XA) | bwamem_pair.c:mem_sam_pe (PE)
 ├─ fastmap.c:main_fastmap ── bwamem_extra.c smem_itr ── bwt.c:bwt_smem1a
 ├─ bwtsw2_main.c ── bwtsw2_aux.c:bsw2_aln ── bwt_lite.c ── bwtsw2_core.c:bsw2_core ── bwtsw2_chain.c ── ksw.c:ksw_extend
 └─ bwtaln.c ── bwtgap.c:bwt_match_gap ── bwt.c (occ)   ;  bwase.c / bwape.c read .sai, load .sa/.pac
```

All mappers share `bwa_idx_load` (`bwa.c:318`), which loads `.bwt`+`.sa` (`bwa.c:271-287`), `.ann/.amb` (`bns_restore`, `bntseq.c:168`) and optionally `.pac` (`bwa.c:307-312`). `bwa aln` is the exception: it loads only `.bwt` (`bwtaln.c:171-175`).

## Entry points

| Command | Entry | Notes |
|---|---|---|
| `bwa index` | `main.c:103` -> `bwa_index` `bwtindex.c:209` -> `bwa_idx_build` `bwtindex.c:255` | options `-a is|bwtsw|rb2`, `-p`, `-b block`, `-6` |
| `bwa mem` | `main.c:111` -> `main_mem` `fastmap.c:141` | getopt string `fastmap.c:159`; pipeline `fastmap.c:395` |
| `bwa fastmap` | `main.c:110` -> `main_fastmap` `fastmap.c:408` | prints SMEMs, loads BWT+BNS only (`fastmap.c:444`) |
| `bwa bwasw` (= `bwtsw2`, `dbwtsw`) | `main.c:107-109` -> `bwa_bwtsw2` `bwtsw2_main.c:11` -> `bsw2_aln` | marked DEPRECATED in usage (`main.c:70`) |
| `bwa aln` | `main.c:104` -> `bwa_aln` `bwtaln.c:230` -> `bwa_aln_core` `bwtaln.c:159` | writes binary `.sai` |
| `bwa samse` / `sampe` | `main.c:105-106` -> `bwase.c:579` / `bwape.c` | `.sai` -> SAM |
| `bwa shm` | `main.c:112` -> `main_shm` `bwashm.c:183` | stage index in POSIX shm; `mem` tries shm first (`fastmap.c:362`) |
| `fa2pac`, `pac2bwt`, `pac2bwtgen`, `bwtupdate`, `bwt2sa` | `main.c:98-102` | manual index pipeline; `bwt2sa -i` sets SA interval |
| C API | `bwamem.h:130-207` (`mem_opt_init`, `mem_process_seqs`, `mem_align1`, `mem_reg2aln`, SMEM iterator) | `example.c` shows use |

## Design decisions the authors made (and why, if stated)

- One FM-index over forward + reverse-complement genome, 64-bit ranks ("works for >4 GB references; twice the memory of 32-bit", `bwa.1:789-793`). Bidirectional extension (`bwt_extend`, `bwt.c:262-275`) exploits this.
- Occ checkpoints are interleaved with the BWT text in a single array (`bwt_bwtupdate_core`, `bwtindex.c:150-172`) so one cache line serves both; the layout is hard-coded for OCC_INTERVAL=128 (`bwt.h:73-75`).
- Seeding by SMEMs (from the fermi assembler paper, `bwa.1:865-870`), then re-seeding and a "LAST-like" third pass (`bwamem.c:169`), chosen for speed on ~0.1-1% error reads.
- Scoring theory in the header comment `bwamem.c:49-70`: score/mismatch ratio encodes an assumed error rate; MAPQ derives from score difference in phred units (6.02 per unit score at A=1).
- Z-drop instead of X-drop "doesn't penalize gaps in one of the sequences" (`bwa.1:181-194`).
- `N` in the reference is replaced by a random base with fixed seed 11 (`bntseq.c:266,295`); positions are recorded in `.amb` so hits spanning N runs can be flagged.
- README (top) now points users to minimap2 (long reads), minibwa (short reads), bwa-mem2 / bwa-mem3 for identical legacy output; this repo is maintenance-only.
