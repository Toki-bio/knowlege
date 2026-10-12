---
name: ncbi-blast-core
description: Internals of the NCBI BLAST+ engine (C core in algo/blast/core plus the C++ api/blastinput layers) for nucleotide search - blastn, megablast, dc-megablast, blastn-short, rmblastn. Use when reading, modifying, debugging or explaining blastn/megablast source, parameters, lookup tables, seeding, X-drop extension, HSP filtering, Karlin-Altschul statistics, -task presets, -subject vs -db, or the megablast database index; especially when diagnosing missed hits for short, diverged, high-copy queries (SINEs, repeats).
---

# NCBI BLAST+ nucleotide engine internals

> Type: reference (facts verified against the pinned commit) with marked `notes`.

## Pinned source

| Field | Value |
|---|---|
| Repo | https://github.com/ncbi/ncbi-cxx-toolkit-public |
| Version / tag | BLAST+ 2.17.0 (`kBlastMajorVersion/Minor/Patch` `src/algo/blast/core/blast_engine.c:81-83`) |
| Commit | `ddcf2de` (full SHA in references/architecture.md), dated 2026-10-10 |
| Digested on | 2026-10-11 |
| Checkout | **sparse, shallow, blob-filtered** clone at `/home/user/src/ncbi-blast`: only `src,include/algo/blast/{core,api,blastinput}`, `src/app/blast`, `src/algo/blast/dbindex`. Unit tests for core, `objtools/blast/seqdb_reader`, `include/algo/blast/dbindex` are NOT present. |
| Language / build | C (core), C++ (api, blastinput, app); CMake (`CMakeLists*.txt`) and NCBI `Makefile.*` |
| License | NCBI Public Domain Notice (`LICENSE`) |

## What it does in three sentences

The BLAST engine indexes the *query* (not the database) into a lookup table of
k-mers, streams every subject sequence from a `BlastSeqSrc`, finds exact word
hits, extends them ungapped with an X-drop, then gapped (greedy or affine DP),
and saves HSPs that beat a raw score cutoff derived from the E-value threshold.
Every behaviour is controlled by option structs (`LookupTableOptions`,
`BlastInitialWordOptions`, `BlastExtensionOptions`, `BlastScoringOptions`,
`BlastHitSavingOptions`) filled by the C++ task presets and `blastn` flags.

## Where to look (read on demand)

| Question | File |
|---|---|
| Directory map, who calls whom, entry points | `references/architecture.md` |
| The structs everything revolves around | `references/data-structures.md` |
| Numbered trace of one blastn run (setup -> prelim -> traceback) | `references/control-flow.md` |
| Lookup tables, word size -> table width/stride, templates, masking, word-frequency filter | `references/lookup-and-seeding.md` |
| Two-hit window, ungapped X-drop, ungapped cutoff, gapped extension, HSP saving, culling/limits | `references/extension-and-hsp-saving.md` |
| reward/penalty tables, allowed gap costs, Karlin-Altschul, search space, E-value, length adjustment | `references/scoring-and-statistics.md` |
| Exact `-task` presets and every default number | `references/task-presets.md` |
| `-db` vs `-subject`, seqsrc, memory cost, megablast dbindex | `references/database-and-index.md` |
| Inputs/outputs and the flags that change them | `references/io-formats.md` |
| Why short diverged repeats are missed; surprising defaults | `references/pitfalls.md` |
| Build, run, debug | `references/build-run-test.md` |

Read only the file that answers the current question. Start with
`architecture.md` when the task is unfamiliar; start with `pitfalls.md` when
the question is "why does blastn miss my hits".

## Ten facts worth knowing before touching the code

1. `blastn` the *program* defaults to `-task megablast` (word 28, 1/-2, greedy) (`src/algo/blast/blastinput/blastn_args.cpp:54-55`). `-task blastn` is word 11, 2/-3, gap 5/2, DP.
2. The lookup table is built from the query; the database is scanned, never indexed, unless `-use_index` (megablast only). Memory is query-bound (`src/algo/blast/core/blast_nalookup.c:46-196`).
3. `-word_size W` means "W consecutive exact matches required"; the table may use a narrower width `lut_word_length` and a scan stride `W - lut_word_length + 1`; hits are then exact-extended to W (`na_ungapped.c:1060-1146`).
4. There is **no** word-frequency / repeat-dropping filter for blastn or megablast. `db_filter` (`-limit_lookup`, max count 30) exists only for magicblast with word size >= 16 (`blast_options.c:1417-1429`, `blast_args.cpp:1495-1502`).
5. Default two-hit window is 0 (one-hit) for blastn and megablast; 40 for dc-megablast (`include/algo/blast/core/blast_options.h:58-61`).
6. In a gapped search the ungapped-extension cutoff is `gap_trigger` = 27 bits converted with the *ungapped* lambda/K (`blast_parameters.c:336-369`), not the E-value; `-xdrop_ungap` 20 bits.
7. Reward/penalty pairs and gap costs are table-driven (`blast_stat.c:611-725`, `3239-3373`); unsupported gap costs below the table max are silently replaced (`3403-3453`), unsupported reward/penalty is an error.
8. For blastn-type programs the E-value is plain Karlin-Altschul on `eff_searchsp = (dblen - N*L)(qlen - L)` (`blast_setup.c:826-839`); the Spouge/Gumbel FSC is not used for blastn (`blast_setup.c:54-62`), and there is no composition-based correction.
9. `-perc_identity`, `-qcov_hsp_perc`, `-max_hsps`, `-subject_besthit` are applied **after** traceback (`blast_traceback.c:606-620, 676-678, 837-867, 1764`); `-culling_limit`/best-hit are HSP-stream writers chosen in `setup_factory.cpp:297-357`.
10. `-subject file.fa` works without `makeblastdb`: subjects are loaded whole into RAM as ncbi2na + ncbi4na (`blast_setup_cxx.cpp:742-744`), and (default, non `BL2SEQ_LEGACY`) are treated as a database for statistics (`blast_app_util.cpp:206-210`, `seqsrc_multiseq.cpp:176-182`).

## Open questions / unverified

- Multi-threaded (`-num_threads`) behaviour with `-subject` (CMultiSeqInfo has a copy ctor, so it should work) not verified by running.
- `BlastHitSavingOptions::min_hit_length` has no blastn CLI flag in this checkout (grep of blastinput found none).
- 8th column of the `blastn_values_*` tables (87..100) is not consumed by any function I read; probably "theta"/percent identity bound. Unused for blastn.
- The dbindex on-disk format (`include/algo/blast/dbindex/*.hpp`) is absent from the sparse checkout; only `README.usage` and `.cpp` were read.
- Nothing here was executed; all numbers under `notes` are arithmetic on the cited constants.

## Change log

- 2026-10-11: initial digest at `ddcf2de` (BLAST+ 2.17.0 era), focus on nucleotide search for SINE consensus vs whole genome.
