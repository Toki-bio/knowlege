# Control flow

## `blastn -query q.fa -db genome` (or `-subject genome.fa`)

1. Parse args: `CBlastnAppArgs` ctor registers argument groups
   (`src/algo/blast/blastinput/blastn_args.cpp:46-121`); default task
   `megablast` (:54-55). `x_CreateOptionsHandleWithTask` ->
   `CBlastOptionsFactory::CreateTask` (`api/blast_options_handle.cpp:337-382`)
   builds the preset, then each `IBlastCmdLineArgs::ExtractAlgorithmOptions`
   overrides it (e.g. `CGenericSearchArgs` `blast_args.cpp:249-323`,
   `CNuclArgs` :673-696, `CFilteringArgs` :386-468, `CWindowSizeArg` :484-498,
   `CHspFilteringArgs` :3331-3346, `CMbIndexArgs` :3375-3411,
   `CFormattingArgs` :2868-2990 which sets `hitlist_size` from
   `-max_target_seqs`/`-num_alignments` :2921-2975).
2. `CBlastnApp::Run` (`src/app/blast/blastn_app.cpp:111-170`): decide
   MT-by-query vs MT-by-DB (`s_GetMT1Cutoffs` :93-109: word_size < 10 -> DB
   <= 150 Mb and query file > 1.5 MB switches to query splitting), then
   `x_RunMTBySplitDB` (:172).
3. `InitializeSubject` (`blast_app_util.cpp:163-225`): `-subject` ->
   `CLocalDbAdapter(subjects, opts, dbscan_mode=true)` unless env
   `BL2SEQ_LEGACY` (:206-210); `-db` -> `CLocalDbAdapter(CSearchDatabase)`
   (:211-220). `MakeSeqSrc` (`api/local_db_adapter.cpp:115-137`) creates the
   `BlastSeqSrc`: SeqDB (`seqsrc_seqdb.cpp`) or in-memory multiseq
   (`seqsrc_multiseq.cpp:558-574`).
4. Query filtering (C++): `CBlastQuerySourceOM::x_CalculateMasks`
   (`api/blast_objmgr_tools.cpp:180-229`) runs dust (default on,
   `Blast_FindDustFilterLoc`, 20/64/1), optional repeats DB
   (`-filtering_db`) and WindowMasker; dust intervals are unioned with
   user/lowercase masks (`api/dust_filter.cpp:92-145`). Lowercase masks come
   from the FASTA reader only when `-lcase_masking`
   (`blast_fasta_input.cpp:397-412`, `cmdline_flags.cpp:115-116`).
5. `SetupQueries` (`api/blast_setup_cxx.cpp:505-660`): blastna encoding,
   both strands as two contexts, masks -> `BlastMaskLoc` stored in
   `seqblk->lcase_mask` (:656-657).
6. `BLAST_MainSetUp` (`core/blast_setup.c:563-662`): compute mask locations
   (:603-608); if `mask_at_hash` (soft masking, default TRUE) the sequence is
   **not** overwritten (:621-623), masks only exclude words from the lookup
   table via `BLAST_ComplementMaskLocations` -> `lookup_segments` (:634-637);
   `BlastSetup_ScoreBlkInit` builds the 16x16 matrix and ungapped K-A blocks.
7. `LookupTableWrapInit_MT` (`core/lookup_wrap.c:61-204`): estimates table
   entries from `lookup_segments` (:127-128), `BlastChooseNaLookupTable`
   (:129-131; `blast_nalookup.c:46-196`), builds MB / small / standard table
   (:133-158). See lookup-and-seeding.md.
8. `BLAST_GapAlignSetUp` (`core/blast_setup.c:888-998`): db length and number
   of sequences from seqsrc (:905-934), `BLAST_CalcEffLengths` (:699-850),
   `BlastScoringParametersNew`, `BlastExtensionParametersNew` (raw X-drops),
   `BlastHitSavingParametersNew` (raw E-value cutoffs), `BLAST_GapAlignStructNew`.
9. HSP writer: `CSetupFactory::CreateHspWriter` (`api/setup_factory.cpp:297-357`)
   -> collector (default), culling (`-culling_limit`, prelim limit +3 :338-340)
   or best-hit (prelim overhang forced to 0 :327).
10. Megablast index: `CSetupFactory::InitializeMegablastDbIndex`
    (`setup_factory.cpp:429-463`) only if `-use_index`; disabled with a
    warning for non-blastn programs, discontiguous templates, or word size
    below index minimum.
11. `BLAST_PreliminarySearchEngine` (`core/blast_engine.c:1335-1605`):
    `BlastInitialWordParametersNew` (:1372-1374) computes ungapped cutoffs;
    `s_BlastSetUpAuxStructures` (:971-1073) picks `BlastNaWordFinder` or
    `MB_IndexedWordFinder` (:1022-1032) and `BLAST_GetGappedScore` (:1072).
    Loop over subjects (:1402-1580): if `db_length == 0` (legacy bl2seq
    mode) recompute eff. lengths per subject (:1430-1439);
    `s_BlastSearchEngineCore` (:702-920).
12. `s_BlastSearchEngineCore` -> `s_BlastSearchEngineOneContext` (:423-603):
    split subject into 5 Mb chunks (`s_GetNextSubjectChunk` :221-312);
    `WordFinder` (:491-499) fills `init_hitlist` (ungapped HSPs sorted by score);
    `GetGappedScore` (:523-527); purge common endpoints (:549); sort (:560);
    merge chunk results (:583-586); subject best-hit if requested (:586-591).
    Back in Core: E-values (`Blast_HSPListGetEvalues` :881-888, or
    `BLAST_LinkHsps` :870-873 when sum stats) and reap by `prelim_evalue`
    (:895-898). `BlastHSPStreamWrite` (:1553) hands the list to the writer.
13. Traceback: `BLAST_ComputeTraceback_MT` (`core/blast_traceback.c:1436-1787`)
    sets `gap_x_dropoff = gap_x_dropoff_final` (:1470, :1424-1433), fetches
    each subject again (uncompressed ncbi4na for ambiguity-aware rescoring)
    and calls `Blast_TracebackFromHSPList` (:259-731) per HSP list; then
    mask-level (rmblastn), `s_FilterBlastResults` (:1764-1766:
    `-max_hsps`, `-qcov_hsp_perc`, `-subject_besthit`), sort by E-value,
    `s_BlastPruneExtraHits` to `hitlist_size` (:1776).
14. Formatting (`src/app/blast/blast_app_util.cpp`, `objtools/align_format`, not in checkout).

## `Blast_TracebackFromHSPList` detail (`blast_traceback.c:259-731`)

1. For each prelim HSP (score order): skip if contained in an already
   traced HSP (`BlastIntervalTreeContainsHSP` :384-386, `min_diag_separation`).
2. Choose start point: `BlastGetStartForGappedAlignmentNucl` (:455-458;
   `blast_gapalign.c:3323`: look for a run of >= 10 identities).
3. Align: greedy traceback (`BLAST_GreedyGappedAlignment` :503-506) for
   megablast, else `BLAST_GappedAlignmentWithTraceback` -> `ALIGN_EX`
   (:507-512, `blast_gapalign.c:4553-4719`).
4. rmblastn complexity adjustment (:521-583).
5. `Blast_HSPUpdateWithTraceback`; identity + length test (non-greedy only
   here :606-620); add to interval tree.
6. After the loop: purge shared endpoints (:650-651), `Blast_HSPReevaluateWithAmbiguitiesGapped`
   (:670-675, rescoring with real ambiguity scores and cutting to best
   sub-segment >= cutoff), `Blast_HSPTestIdentityAndLength` (:676-678),
   blastn second endpoint purge (:682-684), sort, containment purge (:692-704),
   `s_HSPListPostTracebackUpdate` (:199-253): E-values, reap by `-evalue`,
   bit scores.

## Threading / parallelism model

- Prelim: `CBlastPrelimSearch::x_LaunchMultiThreadedSearch`
  (`api/prelim_stage.cpp:142-...`) runs N threads each iterating the shared
  `BlastSeqSrc` chunk iterator (`BlastSeqSrcIteratorNewEx(num_seqs/100)`,
  `blast_engine.c:1401`). Traceback MT in `BLAST_ComputeTraceback_MT`.
- MT-by-query (`blastn_app.cpp:x_RunMTBySplitQuery`) splits the query batch
  instead (`CMTArgs::eSplitByQueries`).

## Error handling conventions

Core returns `Int2` status (0 ok, `BLASTERR_*` codes in `blast_message.h`);
`Blast_MessageWrite` collects warnings/errors; the C++ layer throws
`CBlastException` (`s_BlastMessageToException`). Option validation:
`BLAST_ValidateOptions` (`blast_options.c:1751`) -> `LookupTableOptionsValidate`
(:1284), `BlastScoringOptionsValidate` (:862), etc.
