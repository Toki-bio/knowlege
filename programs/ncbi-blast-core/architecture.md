---
title: "ncbi-blast-core: Architecture"
type: reference
generated: true
---

> Part of the [ncbi-blast-core](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/ncbi-blast-core/references/architecture.md)

{% raw %}
# Architecture

Full commit SHA: `ddcf2de2827fc8ee3eb0d9823718f3f1393db1be`
(https://github.com/ncbi/ncbi-cxx-toolkit-public, commit dated 2026-10-10;
`kBlastMajorVersion=2, Minor=17, Patch=0` at `src/algo/blast/core/blast_engine.c:81-83`).

Sparse checkout: `src/algo/blast/core`, `include/algo/blast/core`,
`src/algo/blast/api`, `include/algo/blast/api`, `src/app/blast`,
`src/algo/blast/blastinput`, `include/algo/blast/blastinput`,
`src/algo/blast/dbindex`. All paths below are relative to the repo root.

## Directory map

| Path | Role | Key files (nucleotide search) |
|---|---|---|
| `src/algo/blast/core/` | C engine, no NCBI toolkit dependency ("blast" lib, `README`) | `blast_engine.c` (prelim loop), `blast_nalookup.c` (query k-mer tables), `blast_nascan.c` (subject scanners), `na_ungapped.c` (word finder, two-hit, ungapped X-drop), `blast_gapalign.c` (greedy + affine DP, traceback), `greedy_align.c`, `blast_hits.c` (HSP lists, E-values, reaping), `blast_traceback.c`, `blast_stat.c` (K-A tables), `blast_parameters.c` (bits -> raw cutoffs), `blast_options.c` (defaults, validation), `blast_setup.c` (masking, eff. lengths), `blast_filter.c` (dust/mask locations), `link_hsps.c` (sum statistics), `hspfilter_{collector,culling,besthit}.c` (HSP writers), `blast_seqsrc.c` (subject source ADT), `lookup_wrap.c` |
| `include/algo/blast/core/` | Headers = the option/parameter structs | `blast_options.h`, `blast_parameters.h`, `blast_nalookup.h`, `blast_hits.h`, `blast_def.h`, `blast_extend.h`, `blast_stat.h`, `blast_seqsrc.h`, `mb_indexed_lookup.h` |
| `src/algo/blast/api/` | C++ wrapper: option handles, task presets, query/subject setup, seqsrc implementations, dbindex glue | `blast_options_handle.cpp` (task factory), `blast_nucl_options.cpp`, `disc_nucl_options.cpp`, `blast_options_cxx.cpp`, `blast_options_local_priv.*`, `blast_setup_cxx.cpp`, `seqsrc_seqdb.cpp`, `seqsrc_multiseq.cpp`, `local_db_adapter.cpp`, `local_blast.cpp`, `prelim_stage.cpp`, `traceback_stage.cpp`, `setup_factory.cpp`, `blast_dbindex.cpp`, `dust_filter.cpp`, `repeats_filter_cxx.cpp`, `winmask_filter.cpp` |
| `src/algo/blast/blastinput/` | Command-line argument classes shared by all BLAST+ apps | `blast_args.cpp` (every flag), `blastn_args.cpp`, `rmblastn_args.cpp`, `cmdline_flags.cpp` (flag names + defaults), `blast_input_aux.cpp` (query batch sizes, `ReadSequencesToBlast`), `blast_fasta_input.cpp` (lowercase masks) |
| `src/app/blast/` | The executables | `blastn_app.cpp`, `blast_app_util.cpp` (`InitializeSubject`) |
| `src/algo/blast/dbindex/` | MegaBLAST database index library + `makembindex` | `README.usage`, `dbindex_factory.cpp`, `dbindex_search.cpp`, `makeindex/mkindex_app.cpp` |

## Module dependency sketch (top-down, blastn)

```
blastn_app.cpp  CBlastnApp::Run                       (src/app/blast/blastn_app.cpp:111)
  -> CBlastnAppArgs (blastinput)   parses flags -> CBlastOptionsHandle (task preset + overrides)
  -> InitializeSubject (blast_app_util.cpp:163)  -db -> CLocalDbAdapter(CSearchDatabase)
                                                  -subject -> CLocalDbAdapter(IQueryFactory, dbscan=true)
  -> CLocalBlast::Run (api/local_blast.cpp:167)
       -> CBlastPrelimSearch (api/prelim_stage.cpp)
            -> SetupQueries (api/blast_setup_cxx.cpp)  masks: dust/lcase/windowmasker -> BlastMaskLoc
            -> BLAST_MainSetUp        (core/blast_setup.c:563)  lookup_segments, score block
            -> LookupTableWrapInit    (core/lookup_wrap.c:61)   query k-mer table
            -> BLAST_GapAlignSetUp    (core/blast_setup.c:888)  eff. lengths, cutoffs, gap_align
            -> CSetupFactory::CreateHspWriter (api/setup_factory.cpp:297) collector|culling|besthit
            -> BLAST_PreliminarySearchEngine (core/blast_engine.c:1335)
                 for each subject: s_BlastSearchEngineCore (:702) -> s_BlastSearchEngineOneContext (:423)
                    -> BlastNaWordFinder (core/na_ungapped.c:1605)  scan + exact extend + ungapped X-drop
                    -> BLAST_GetGappedScore (core/blast_gapalign.c:3664) greedy | DP score-only
                    -> Blast_HSPListGetEvalues / reap (core/blast_hits.c:1811, blast_engine.c:643)
                 -> BlastHSPStreamWrite -> writer (hspfilter_*.c)
       -> CBlastTracebackSearch -> BLAST_ComputeTraceback_MT (core/blast_traceback.c:1436)
            -> Blast_TracebackFromHSPList (:259) per subject: traceback alignment, identity/coverage
               filters, containment purge, final E-values
```

## Entry points

| Command / mode | Entry | Notes |
|---|---|---|
| `blastn` | `CBlastnApp::Run` `src/app/blast/blastn_app.cpp:111-170` | picks MT-by-query vs MT-by-DB (`s_GetMT1Cutoffs` :93-109), then `x_RunMTBySplitDB` :172 |
| `blastn -task ...` | `CBlastOptionsFactory::CreateTask` `src/algo/blast/api/blast_options_handle.cpp:337-382` | see task-presets.md |
| `rmblastn` | `CRMBlastnAppArgs` `src/algo/blast/blastinput/rmblastn_args.cpp:40-121` | adds `CRMBlastNArg` (`blast_args.cpp:526-571`) |
| `makembindex` | `src/algo/blast/dbindex/makeindex/mkindex_app.cpp:66-200` | index for `-use_index` |
| Core, no C++ | `Blast_RunFullSearch` `src/algo/blast/core/blast_engine.c:1727` | prelim + traceback in one call |

## Design decisions the authors made (and why, if stated)

- Query is indexed, subject is streamed in 2-bit (ncbi2na) packed form; scanners
  read 4 bases per byte (`COMPRESSION_RATIO 4`, `blast_def.h:83`).
- Lookup-table width is a cache trade-off chosen per query size
  (`blast_nalookup.c:51-58` comment).
- Long subjects are processed in chunks of `MAX_DBSEQ_LEN` = 5,000,000
  (`blast_gapalign.h:54`) overlapping by `DBSEQ_CHUNK_OVERLAP` = 100
  (`blast_hits.h:192`); HSPs crossing chunk boundaries are merged
  (`blast_hits.c:2857`). Hence `-word_size` must be <= 100 (`blast_options.c:1325-1333`).
- Ungapped nucleotide extension is first done approximately on 4-base
  packed bytes and redone exactly only above 0.8x cutoff
  (`na_ungapped.c:262-357`, `blast_parameters.c:410-414`).
- Statistics for blastn remain classic Karlin-Altschul; Gumbel/Spouge FSC is
  explicitly "not supported for nucleotide search yet" (`blast_setup.c:54-62`).
- rmblastn (RepeatMasker) support is bolted on via `matrix_only_scoring` and
  `complexity_adjusted_scoring` flags, marked `-RMH-` throughout.
{% endraw %}
