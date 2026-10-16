# Pitfalls and gotchas

HitSV v2.0.4, commit `93880367ef76f6440c918f7dd9744dda2659bcdb`. Paths are relative to `src/`
unless they start with a repo-root file name. "Tested" means reproduced on a binary built from
the pinned commit (see `build-run-test.md`).

## Read these first

| # | Gotcha | Evidence | Workaround |
|---|---|---|---|
| 1 | Wrapper script SRS/Hybrid modes call a non-existent `srs_bam_stat`; `stat.json` ends up empty and exit code is 0 | `HitSV_multy_process.sh:145`; dispatch `cpp_lib/get_option_cpp.hpp:367`; registered name `main.cpp:106`. Tested | Use `ngs_bam_stat`, or patch the script |
| 2 | `-I` repeat file is written in FASTA record order but queried by BAM `tid` | `SVcalling_core/fa_stat.hpp:97`; `SV_core_SRS_ins_del.cpp:664` | Build `-I` from the exact FASTA used for alignment; keep BAM header order = FASTA order |
| 3 | Human Alu force-calling is always on; synthetic ~320 bp IMPRECISE INS with position-derived length | `sv_call_main.hpp:175`; `SV_core_SRS_inv_tra.cpp:822`; `SV_core_SRS_ins_del.cpp:1801` | The internal `force_ALU` flag (`NovaSVRst.hpp:1751-1752`) is not written to the VCF; in SRS mode drop INS with `IMPRECISE` and SVLEN 313–327 (`NovaSVRst.hpp:675-676`) |
| 4 | Hybrid "polishing" (`-k`, v2.0.4) edits a temp buffer that is then rebuilt from the unpolished contig: no effect on VCF | `SV_core_SRS_ins_del.cpp:3206-3210` vs `:3357,3366` | Do not expect `-k` to change calls (inferred from code; not run-tested) |
| 5 | Every written LRS record is PASS; 0/0 is never written; use QUAL (20–500 heuristic) | `NovaSVRst.hpp:624,660`; `SV_core_LRS.cpp:2292-2312` | Filter on QUAL / support fields, not FILTER |
| 6 | Small variants (SNP/INDEL) are mixed into the SV VCF by default | `sv_call_main.hpp:135`; demo run output | `-v 0` or filter on `SVTYPE` |
| 7 | "Phased" `0\|1`/`1\|0` is pseudo-random phasing, on by default | `sv_call_main.hpp:133`; `NovaSVRst.hpp:520-537` | `-f 0` unless you run the cohort CSV tools (which need it) |
| 8 | Default `-P asm5` contig realignment is tuned for human | `sv_call_main.hpp:150-153`; `Contig_aligner.hpp:801-813` | `-P asm10`/`asm20` for diverged references |
| 9 | Collapsed repeats/satellites are skipped (depth gates) | `SV_core_SRS_ins_del.cpp:907-914`; `SV_core_LRS.cpp:2320-2321` | None in code; expect blind spots in satellite-rich blocks |
| 10 | Many tools cap at 24 reference sequences (chr1-22,X,Y) | `SVcalling_core/ANALYSIS/analysis.hpp:2857`; `SV_core.hpp:876-882` | Put target scaffolds first in FASTA and BAM header |
| 11 | `ngs_merging` / `ngs_fc` read hard-coded developer paths and do not produce usable output | `NGS_merging.hpp:184-185`; `forceCalling.hpp:873,979` | Do not use; merge with SURVIVOR/Jasmine as the paper did |
| 12 | Documented cohort workflow commands fail as written (`fc_joint_ana`, missing `-B PURE_STR`, MEI-TR script args) | see section D below. Tested | Use the corrected commands in `population-csv-mei-tr.md` |

## Defaults that surprise

| Parameter | Default | Why it matters |
|---|---|---|
| `-m MIN_SV_len` | 25 (`sv_call_main.hpp:134`) | Paper benchmarks use ≥50 bp; filter downstream |
| `-p LRS_preset` | `ONT_Q20` (`sv_call_main.hpp:144`) | Set `HIFI` for PacBio; `ERR_PRONE` disables small-var output (`:218-220`) |
| `-e MAX_read_ERR_RATE` | `"0.04"`, ignored | Dead option (`sv_call_main.hpp:145-146`); real thresholds hard-coded (`SV_core_LRS.cpp:2519-2521`) |
| `-E ed_chr_ID` | 50000000 (`sv_call_main.hpp:126`) | All contigs; `-S/-E` are 0-based BAM header indices, not names |
| Block size | 5 Mb, no overlap (`RefHandler.hpp:30`) | Comment says 2 Mb; FC_BED regions crossing a block edge are dropped (`SV_core.hpp:889`) |
| Sample name | always `Sample` (`NovaSVRst.hpp:1696`) | Rename before merging cohorts |


## Pitfalls (A: data structures, reference, stats, I/O)

All paths are relative to `src/`. Commit 93880367.

### Reference and stat files

- **The `-I` stat file is keyed by FASTA record index, but queried by BAM `tid`.** `srs_fa_stat` writes `i` = `faidx_iseq` order (`SVcalling_core/fa_stat.hpp:97`). `RefLocalRepeat::search` is called with `region.chr_ID`, the BAM tid (`SV_core_SRS_ins_del.cpp:664`). If the FASTA order differs from the BAM header order (for example after `samtools reheader` or a sorted dict), the repeat tags land on the wrong chromosomes, and nothing warns.
- **The stat file is silently partial.** Only rows with R ≥ 40 are kept (`RefHandler.hpp:123`). Chromosomes with no rows give no hits and no warning (`:136-139`).
- **Exactly 3 tab-separated columns are required** on every line (`xassert`, `RefHandler.hpp:111`). The BED-style output (`fa_stat.hpp:99`) has 4 columns and would abort, but it cannot be produced anyway: `run_ana_bed` is not registered in `main.cpp:105-115`.
- **Soft-masked FASTA changes R.** k-mers are compared as raw strings, case-sensitive. `N_count` only checks for `'N'` (`fa_stat.hpp:17-24`), so lowercase `n` runs are counted as sequence.
- **The window and step must stay 300-aligned.** `search` probes integer positions and expects window starts on multiples of `region_step_size` (`RefHandler.hpp:132-139`), and the file is always written with step 300 (`fa_stat.hpp:35, 49`). The parameters are hard-coded, with no CLI options.
- **A dead fallback for a missing `-I`.** `run` has an "Genome Context information is empty" branch (`sv_call_main.hpp:267-271`). It is unreachable, because `get_option` already rejects `-n` without `-I` (`:197-200`). Note: `-i` in README.md:208 is the *wrapper script's* flag (`HitSV_multy_process.sh:81`), which passes it on as `-I` (`HitSV_multy_process.sh:173`); `call` itself only accepts `-I` (`sv_call_main.hpp:116`).
- **`srs_fa_stat` is single-threaded and uses `std::set<std::string>`.** Each window needs up to 9 passes over about 580 heap-allocated strings once any 20-mer repeats, because `withloop` is never reset (`fa_stat.hpp:68-91`). Expect long runtimes on multi-Gb genomes (inferred, not timed).

### BAM statistics

- **The bundled wrapper script is broken for SRS/Hybrid.** `HitSV_multy_process.sh:145` calls `srs_bam_stat`, which is not a registered command (exact `strcmp` dispatch, `cpp_lib/get_option_cpp.hpp:367`). The binary prints `wrong command`, exits 0 and leaves `stat.json` empty; every `call` then gets `-T` pointing at an empty file (`HitSV_multy_process.sh:173,177`). Tested on the built binary. Fix: replace with `ngs_bam_stat`.
- **Help text names the wrong command.** It says `srs_bam_stat` (`main.cpp:107`; also in `-T` help, `sv_call_main.hpp:158`), but the registered name is `ngs_bam_stat` (`main.cpp:106`).
- **`ngs_bam_stat` does not validate its arguments.** It uses `argv[optind]` and `argv[optind+1]` (`main.cpp:95`), so missing arguments would dereference NULL (inferred).
- **Without `-T`, every `call` invocation recomputes the BAM stats** by sampling the BAM (`sv_call_main.hpp:276-277`, `bam_stat.hpp:140-157`). This happens once per region job in parallel `-S/-E` runs.
- **The read-length estimate is capped.** Reads ≥ 1000 bp are not counted in the histogram or the sum, but they do increase `total_read_number` (`bam_stat.hpp:94-99`). Mixed or long SRS lengths therefore bias `ave_read_len` low.
- **`getBreakPoint_Distribution` assumes `maxInsertLen > 2*analysis_read_length`** (`bam_stat.hpp:166-173`). Short-insert libraries, where the insert is shorter than twice the read length, give a non-positive vector size (inferred).

### Read filters

- **`-e MAX_read_ERR_RATE` is a dead option.** It is parsed into `MAX_read_error_rate` (`sv_call_main.hpp:146`; comment `//todo: not work right now` at `:145`) and never read anywhere. The LRS thresholds are hard-coded at 0.04 / 0.16 (`SV_core_LRS.cpp:2519-2521`).
- **`-e` may also be stored with the wrong type.** Its default is the string `"0.04"`, which selects the `const char*` overload (`get_option_cpp.hpp:171`), but the target field is a `double` (`sv_call_main.hpp:47`). The setter would then write a `char*` into the double field (`get_option_cpp.hpp:107, 117`). This is harmless only because the value is never used.
- **LRS reads are not de-duplicated.** There is no `bam_is_duplicate` test in `SV_core_LRS.cpp` (grep: 0 hits), while SRS reads are de-duplicated (`SV_core_SRS_ins_del.cpp:437`).
- **The SRS MAPQ rescue depends on the `MQ` tag.** It returns 0 when the tag is absent (`clib/bam_file.c:565`). Aligners that do not write `MQ` lose the rescue of MAPQ<20 reads whose mate maps well (`SV_core_SRS_ins_del.cpp:21-23`).
- **The SRS depth estimate counts reads that are later dropped**, such as secondary and unmapped reads (`SV_core_SRS_ins_del.cpp:439-447`, `:476`). This inflates the depth and so the depth-scaled cut-offs (`SV_core.hpp:70-81`).

### Blocks and regions

- **`-S/-E` are BAM header indices, not names.** The help text says "-S 0 -E 23 … human whole-genome chr1 to chrY" (`sv_call_main.hpp:123`). This is a hard-coded human assumption: for other genomes the index range must be looked up in the BAM header. The default `-E 50000000` (`:126`) covers all contigs.
- **The comments state the block size wrongly.** The block is 5 Mb (`RefHandler.hpp:30`), while the comment says "2M bases per block" (`sv_call_main.hpp:298`). Unused leftovers `REF_BLOCK_LEN/NUM`, `HASH_LEN`, `KMER_LEN` and `REF_index` remain (`RefHandler.hpp:26-68`).
- **Contig index limit of 65536.** `RefRegion` stores `chr_ID` through `uint16_t` parameters (`cpp_lib/RefRegion.hpp:25, 30, 39`). Highly fragmented assemblies with more than 65,535 contigs would wrap the IDs.
- **Many tiny contigs are slow.** Each contig is at least one 5 Mb block iteration, and every pipeline call re-runs a BAM region query and handler setup (`SV_core_LRS.cpp:1964-1968`; `SV_core_SRS_ins_del.cpp:420-423`). This is inferred overhead for scaffold-level assemblies.
- **Somatic BND calls abort past 500 contigs.** `xassert(main_BP_ref_chrID < 500)` (`SV_core_SRS_ins_del.cpp:3582`) fires on any genome with more than 500 contigs.
- **The whole chromosome is held in memory.** It is stored as one DNA5 code per byte (`RefHandler.hpp:314-321`). In addition, each `BAM_handler` allocates a fixed 300 MB `storeReadBuff` (`ReadHandler.hpp:609`), so hybrid runs use at least 600 MB of baseline buffers.

### VCF output

- **The `lt50bp` filter is declared but never produced.** Every call with |SVLEN| < 50 is written as `PASS` (`NovaSVRst.hpp:660`, header at `:1673`). The demo `demo/lrs.vcf` shows `lt50bp`, but it comes from an older build (`##source=GC_SV2.1.0`).
- **`LOW_DEPTH` and `Duplication_sv` cannot appear.** Records with GT 0/0, 0|0 or duplication are returned before output (`NovaSVRst.hpp:624`), even though FILTER logic for them exists (`:661-662`).
- **Random phasing is on by default** (`-f` default 1, `sv_call_main.hpp:133`), and it applies to SRS calls too (`SV_core.hpp:2617-2621`). It is not random: all single-haplotype hets from a block get the same `0|1`/`1|0` orientation, because `global_SV_ID` is only advanced during writing (`NovaSVRst.hpp:632-634`; `SV_core.hpp:2554-2564`). Use `-f 0` for unphased output.
- **Small variants are on by default** (`-v 1`, `sv_call_main.hpp:135`). SNP and INDEL records then share the VCF with SVs (`SV_core_SRS_ins_del.cpp:4019-4073`), and they have FILTER=PASS because |SVLEN| < 50. Downstream SV benchmarking must filter `SVTYPE` or use `-v 0`.
- **`##fileDate` is wrong.** It is `1990+tm_year` without zero padding (`NovaSVRst.hpp:1645`), giving e.g. "2115116". `##source` reports `PACKAGE_VERSION "2.0.0"` (`clib/desc.h:2`), not 2.0.4.
- **Deletions over 50 kb are dropped.** SRS genotyping forces 0/0 for them (`NovaSVRst.hpp:1622-1624`), which also means they are never written. So the symbolic `<DEL>` path (`:643`) can only be reached by LRS calls (inferred).
- **No DUP SVTYPE is emitted.** DUP DR signals exist (`sve.hpp:23`), but no `add_to_vector(... "DUP")` call exists. TRA appears only as `SVTYPE=BND;BND_TYPE=TRA`, and the header-declared ALT `TRA`/`INV` symbols are not used as ALT alleles.
- **`QUAL` is a capped heuristic.** It is clamped to [20, 500], forced to 20 for imprecise length, and divided by 10 when U/A > 0.1 (`SV_core_LRS.cpp:2296-2304`). With A = 0 the division `sup_both/sup_alt` is a float divide by zero (`:2298`), which gives inf or NaN; the result is still bounded by the clamps.
- **`GT_cal.hpp` (`GT_handler`) is unused** (only a demo struct in the same file, `GT_cal.hpp:243-247`), and `GT_cal.cpp` is empty. The real genotyping is the simple ratio rules at `SV_core_LRS.cpp:2279-2290` (LRS) and `NovaSVRst.hpp:1507-1528` (SRS).
- **The sample column is always named "Sample"** (`NovaSVRst.hpp:1696`). It is not taken from the BAM `@RG SM`.
- **`fclose(vcf_w)` also runs when the output is stdout** (`SV_core.hpp:2351`). `RefHandler::destory` and `delete sve_h` are never called (`sv_call_main.hpp:311`).
- **`-B PURE_STR` only works together with `-b` FC_BED.** The basic path always passes `output_pure_contig=false` (`SV_core.hpp:2506`). An unknown `-B` value silently becomes `OUT_MODE_UNSET` (`sv_call_main.hpp:237-242`), which skips the k-mer index build (`:304-305`), but VCF records are still written in the basic path.
- **The ALU force-calling sequence is a hard-coded human Alu consensus** (`SV_core_SRS_inv_tra.cpp:822`). It is always on (`force_calling_ALU = true`, `sv_call_main.hpp:175`) and runs whenever SRS assembly yields no call (`SV_core_SRS_ins_del.cpp:2119`). For non-primate genomes this can produce spurious "INS" calls with `IMPRECISE` (inferred risk), and there is no CLI switch to turn it off.

## Pitfalls B: LRS, assembly, realignment, polishing, genotyping, somatic

HitSV v2.0.4, commit 93880367ef76f6440c918f7dd9744dda2659bcdb. Paths are relative to `src/`.

### Bugs and dead code
- **SRS polishing has no effect on calls.** `contig_polsihing_handler.run(bin_contig, …)` edits the member `bin_contig` (`SVcalling_core/SV_core_SRS_ins_del.cpp:3206-3210`). The result is never written back to `target_path_idx_l[...].contig`. Downstream code rebuilds `bin_contig` from the unpolished string: germline at `:3366`, somatic `bin_tc` at `:3650`. As a result, `-k` and the hybrid polishing step in v2.0.4 change nothing in the VCF; they only cost runtime.
- **`GT_cal.hpp` (Bayesian `GT_handler`) is never used.** No reference exists outside the header, and `GT_cal.cpp` is empty (`GT_cal.cpp:1-10`). LRS genotypes come from a 6:1 ratio rule (`SV_core_LRS.cpp:2279-2290`).
- **`-e MAX_read_ERR_RATE` is a no-op.** It is parsed (`sv_call_main.hpp:146`, comment "todo: not work right now" at :145). Hardcoded thresholds apply instead: 0.04, or 0.16 for ERR_PRONE (`SV_core_LRS.cpp:2519-2521`).
- **Somatic size check uses a stale buffer.** The single-support >10 kb test uses `bin_contig.size()`, the member left over from earlier work, not the tumor contig `bin_tc` (`SV_core_SRS_ins_del.cpp:3674-3675`).
- **Somatic: normal assembly failure drops the tumor signal.** If the normal BAM is given and its assembly fails, the `continue` at `SV_core.hpp:2417` skips the tumor event entirely. The event is not reported as TUMOR_UNIQ/UNKNOWN.
- **Initial ksw z-drop reads an unassigned member.** `set_option(..., gap_open2_D + 600, ...)` is evaluated before `gap_open2_D` is set (`Contig_aligner.hpp:801-813`). Because of the memset in the constructor (`SV_core.hpp:1943`), the value is 600. This is harmless only because `setZdrop` overrides it later (`SV_core_SRS_ins_del.cpp:2921`).
- **Long-mode alignment widens the band permanently.** `align_non_splice_super_long` sets band and z-drop to 30000 on the shared aligner (`Contig_aligner.hpp:205`). Normal-mode calls keep that value until the next `setZdrop` (`SV_core_SRS_ins_del.cpp:2921`).
- **Comment/code mismatches.**
  - min LRS read slice is 800, not "1800" (`SV_core_SRS_ins_del.cpp:2961`);
  - the homopolymer k-mer is 40, not "15" (`:3172`);
  - blocks are 5 Mb, not "2M" (`RefHandler.hpp:30` vs `sv_call_main.hpp:298`).

### Paper claims vs code (LRS mode)
- **"indel ≥30 bp"**: the code threshold is `MIN_sv_len`, default **25** via `-m` (`sv_call_main.hpp:134`; `SV_core_LRS.cpp:2067, 2083`). This does not match the paper unless the paper ran with `-m 30`.
- **"clip ≥300 bp"**: the soft clip must be **>1000 bp** to become a signal (`SV_core_LRS.cpp:2042`). For assembly, a read counts as clipped when soft clip > 0 and the read is not full-span (`ReadHandler.hpp:142-152`). Reads clipped >200 bp at both ends are excluded (`SV_core_SRS_ins_del.cpp:1013`). No 300 bp clip threshold exists in the LRS path. The only 300 values are the clustering distance (`SV_core_LRS.cpp:2192`) and the hybrid dedup window (`:2616`).
- **"10 events in 200 bp"**: this matches. The trigger is >10 events **or >35 weighted bases** within a 200 bp window (`SV_core_LRS.cpp:2539, 2547`). It is active only for ONT_Q20 and HIFI (`:1954-1955`).
- **"depth fraction 10%"**: this matches. The support gate is `0.1×block mean depth` (`SV_core_LRS.cpp:2189, 2316`). The depth is a 5 Mb block average (`:2172-2173`), not local.
- **"±500 bp extension"**: this matches for both reference and reads (`SV_core_SRS_ins_del.cpp:2864-2865`). In hybrid mode, an extra SRS read-length is added on the left only (`:2863, 2872`).
- **Local repetitiveness R**: R is not used for LRS detection or assembly. `rlr` is NULL without `-n` (`sv_call_main.hpp:256, 263-271`). Its only LRS use is skipping SRS re-GT in ERR_PRONE at depth <10 (`SV_core.hpp:2529-2536`).
- **SA/split alignments**: in germline mode they are parsed but never used as signals. They feed only the somatic BND builder (`SV_core_SRS_ins_del.cpp:3758`). Germline translocations and inversions from LRS are found only indirectly, through clip signals and the long-range reference-index path, which handles 1-50 kb (`SV_core.hpp:1606-1607`).

### Non-human genomes and satellite/tandem arrays
- **Diploid is hardcoded.**
  - At most 2 final haplotype contigs per region (`SV_core_LRS.cpp:1265-1269`).
  - Genotypes are limited to 0/0, 0/1, 1/1 (`:2279-2290`).
  - Same-site 0/1+0/1 on two haplotypes is merged to 1/1 (`:2444-2465`).
  - Polyploids and pools are not modeled.
- **ASM preset keeps one sequence per region** (`SV_core_SRS_ins_del.cpp:1253-1274`). A diploid assembly given as one BAM cannot yield 0/1 calls from both haplotypes in one region. Map each haplotype separately.
- **Repeats are invisible to the reference anchor index.**
  - The index is 16-mers sampled every 3 bp, unique only (`RefHandler.hpp:209, 283, 297`). Long-range (clip) contigs in satellite arrays get no clusters (needs >60 hits, `SV_core_LRS.cpp:1819`) and are dropped as "No_sv_can_be handler" (`:1882-1884`).
  - SDP ignores unitigs repeated >30× in a path (`:630`).
- **Expansion/contraction size is limited.**
  - Reads are sliced to the window plus clips, and soft clips ≥100 kb are not loaded (`ReadHandler.hpp:810`).
  - Contigs come only from reads spanning the ±500 bp-padded region; contig selection requires ≥1 full-cover read (`SV_core_LRS.cpp:1220-1222`). An array change longer than the reads cannot be assembled end-to-end.
  - Long-range resolution handles 1-50 kb (`SV_core.hpp:1606-1607`).
- **Single-read evidence in repeats is filtered.**
  - support==1 calls whose ALT is >40% low-complexity 10-mers are dropped unless the 10-mers recur in the flanks (`SV_core_SRS_ins_del.cpp:4142-4241`).
  - Single-support clip contigs >10 kb are skipped (`:3392-3395`).
- **Fragmented calls.** Adjacent CIGAR gaps in a contig alignment are not summed, because the summing code is disabled (`Contig_aligner.hpp:415, 434, 454`; the TGS extractor at `:495-555` has none). A tandem-repeat change can be emitted as several INS/DEL records. Run `bcftools norm`/`truvari collapse` downstream.
- **Depth gates use a 5 Mb block mean** (`SV_core_LRS.cpp:2172-2173`). In heterochromatic or satellite blocks with mapping dropouts, `0.1×mean` can miss or admit regions. Regions with more loaded reads than `MIN(MAX(5×depth,100),1000)` are skipped (`SV_core_LRS.cpp:2320-2321`; `SV_core_SRS_ins_del.cpp:2969-2973`). Collapsed satellites with excess depth are therefore silently skipped.
- **The k-mer DBG has no reference k-mers in LRS mode** (`DBGAssembler.hpp:3178-3179`). With k=20 (10 for ERR_PRONE, `SV_core_SRS_ins_del.cpp:3001-3003`), satellite monomers shorter than about 20 bp collapse into cycles. Cycles are flagged `isRepeat` but the flag is ignored in the LRS path; resolution relies on per-read threading.
- **Memory.** Each chromosome is loaded whole as 1 byte/base (`RefHandler.hpp:314-321`). Somatic tumor assembly can take `MAX(10×depth,100)` reads (`SV_core.hpp:2378`). The VCF line buffer is a fixed 1 MB filled with `sprintf` (`SV_core.hpp:2312-2313`; `NovaSVRst.hpp:628-669`); an INS ALT plus contig/CIGAR INFO near 1 MB would overflow it (unverified at runtime).
- **Every emitted LRS record is PASS.**
  - 0/0 records are not written (`NovaSVRst.hpp:624`);
  - |SVLEN|<50 is always PASS (`:660`);
  - DR/VNTR filters default to pass (`:385`).
  - Filter on QUAL (range 20-500, `NovaSVRst.hpp:565-566`), not FILTER.
- **Random phasing is on by default** (`-f 1`, `sv_call_main.hpp:133`). It is called unconditionally when `random_phasing` is true (`SV_core.hpp:2554-2557`) and turns 0/1 into 0|1 or 1|0. These pipes are not real phasing.

## Pitfalls C: SRS / hybrid / ngs_merging / ngs_fc

HitSV v2.0.4, commit 93880367. Paths are relative to `src/SVcalling_core/` unless noted.

### Human-specific or biasing for non-human genomes

- **A human Alu consensus is compiled in, and force-ALU calling is always on.** It has no CLI switch (sv_call_main.hpp:175; consensus at SV_core_SRS_inv_tra.cpp:822).
  - For a reptile or other mammal genome it will rarely fire, because it requires an exact anchor k-mer match (inv_tra.cpp:901-911).
  - When it does fire, it emits a synthetic ~320 bp IMPRECISE INS (SV_core_SRS_ins_del.cpp:1801, 1960-1965).
  - It is not adapted to lineage-specific SINEs. No other MEI consensus (L1, SVA, squamate SINEs) exists.
- **Synthetic Alu insertions are not real sequence.** The ALT length is `320 + (pos%15 - 7)`, which is position-derived pseudo-random (ins_del.cpp:1801). The ALT string is taken from the contig padded with N (ins_del.cpp:1810-1857). Do not use force-ALU INS sequences or lengths for TE analysis. The `force_ALU` flag is internal (NovaSVRst.hpp:1751-1752) and is not written to the VCF; such records appear as `SVTYPE=INS` with `IMPRECISE` (NovaSVRst.hpp:675-676) and SVLEN 313–327 (320 + pos%15 − 7).
- **The Alu realignment test is length-hardwired**: the CIGAR must contain an I of 280-350 bp (Contig_aligner.hpp:697).
- **Poly-A/T heuristics assume Alu/L1-like tails.**
  - A clipped read is discarded as a signal read when the 20 *aligned* bases next to either clip are >18 A or >18 T, i.e. the clip sits on a reference poly-A/T run (ReadHandler.cpp:23-44, called at :83). The clipped bases themselves are not inspected.
  - Insert reconstruction joins contigs through a >=15 bp homopolymer (SV_core_SRS_combine_contigs.cpp:54, 245).
  - SINEs or other insertions without long A/T tails get no tail-joining help. Insertions with tails lose the clip-signal reads that end in the tail.
- **The contig realignment default is `asm5`, labelled "Recommended for Human".** Use `-P asm10`/`asm20` for diverged species (sv_call_main.hpp:150-153).
- **chrY / 24-contig assumptions** appear only in tools and help text:
  - help example `-S 0 -E 23` (sv_call_main.hpp:123)
  - `ngs_merging` caps at 24 contigs (NGS_merging.hpp:171, 239)
  - The `call` path itself iterates all BAM header targets (RefHandler.hpp:243-258).
- **No hardcoded insert size or read length in `call`.** They come from `-T` JSON or on-the-fly stats (bam_stat.hpp:122-138; sv_call_main.hpp:273-277).
  - Read-length detection ignores reads >=1000 bp (bam_stat.hpp:75, 98).
  - It picks the length carried by >60% of reads (bam_stat.hpp:110-113), else the mean (117-118).
  - `ngs_fc` hardcodes isize 1000 and readlen 150 (forceCalling.hpp:709-710, 889).

### Correctness and behaviour gotchas

- **The `-I` file is indexed by FASTA order, but looked up by BAM tid.** fa_stat writes `i` (FASTA index) as chrID (fa_stat.hpp:52, 97). `RefLocalRepeat::search` is queried with BAM tid (ins_del.cpp:662-664; SV_core.hpp:2531). If the BAM header contig order differs from the FASTA (e.g. a sorted or renamed header, or a subset BAM), repeat windows land on the wrong chromosomes, silently.
- **The R threshold is `>=40`, not `>40`** (RefHandler.hpp:123). R only takes values 20, 40, ..., 180 (fa_stat.hpp:36-38, 70), so R=20 windows are dropped and all others become VNTR candidates (ins_del.cpp:663-671). Repeat-rich genomes (e.g. satellite-dense squamate chromosomes) will generate very many candidates. Each is a dBG run with up to 60000 reads (ins_del.cpp:2238), so runtime grows sharply.
- **The `fa_stat` repeat scan never resets `withloop`** across k values (fa_stat.hpp:68-92). R reports the k at which a repeat was last seen before the loop could not break. Treat R as a coarse proxy, not exact maximal repeat length.
- **The paper's "max(10% depth, 2)" read threshold is really `max(3, 0.10*depth)`**, plus `max(2, 0.08*depth)` on score (SV_core.hpp:72-75). Clusters seen on one side only need 2x both (ins_del.cpp:586). There is an earlier gate of >=3 signals and score > clamp(0.4*depth,16,29) (SV_core.hpp:77-79; ins_del.cpp:315-317). Low-coverage SRS (<=10x) will lose many calls.
- **Depth is block-local**: read count * readlen / 5 Mb block (ins_del.cpp:476; SV_core.hpp:2591). It is not the JSON depth. Blocks that are mostly N or unplaced scaffolds deflate depth and loosen thresholds. Short scaffolds count as one block each (RefHandler.hpp:245-250).
- **High-depth skip.** Assembly is skipped when local mean depth is >5x or a 32-bp bin is >15x the genome average (ins_del.cpp:907-914; normalisation SV_core.hpp:158). Collapsed repeats and satellites in the reference are therefore never assembled.
- **Duplications are never called from SRS.** DR/DUP signals are clustered (ins_del.cpp:550-552, 626-639) but no caller reads `SRS_sve[DR][DUP]`.
- **SRS TRA is pair counting only.** `SRS_assembly_and_genetyping_BND` does no assembly (inv_tra.cpp:76-158). It is stored as "polished" (inv_tra.cpp:202), so the IMPRECISE filter at SV_core_SRS.cpp:459 never drops anything. The assembly function `SRS_assembly_and_get_breakpoints_TRA` (inv_tra.cpp:360) is never called.
- **SRS INV and DEL are capped at 50 kb** (SV_core_SRS.cpp:287, 378, 496-497). DR pairs with |isize| >= 50000 on the same chromosome get no SV type (ins_del.cpp:545; SV_core.hpp:67). The INV/TRA branches handle orientation and chromosome tests only.
- **`-R` (realigned reads) is dead**: `if(false) SRS_load_realignment_read_from_file()` (SV_core.hpp:2590).
- **`MIN_TL_read_number` gate is dead**: computed but never compared meaningfully (SV_core_SRS.cpp:523-529).
- **No `-L` means a slow fallback, not a failure.** It rescans the main BAM per mate cluster (ReadHandler.hpp:388-389, 466-530).
- **TL reads are not used for TRA calling**, only for insert reconstruction (ins_del.cpp:1402; combine_contigs.cpp:569-721).
- **Short-insert libraries can break the breakpoint model.** `max_probability_size = maxInsertLen - 2*readlen` is passed to `resize` (bam_stat.hpp:166-168), and `insert_region_len` uses the same formula (SV_core.hpp:65). If the 99th-percentile insert is <= 2x read length (overlapping 2x150 or 2x250 libraries), this goes to <=0. Behaviour there is unverified, with a likely crash or a degenerate DR model. Check the `-T` JSON before running.
- **The DBG "max visit 9" in the paper is tier 2.** The first try uses 13, then 9, 5, 1 if too complex (DBGAssembler.hpp:3345). The haplotype-pair "likelihood" is a depth-consistency heuristic (DBGAssembler.hpp:1556, 1634-1650), not a probability.
- **The VNTR caller discards a site's calls if it yields both INS and DEL** (ins_del.cpp:2314-2332). Heterozygous expansion/contraction pairs at the same locus are lost.

### Hybrid-mode gotchas

- **SRS calls within 300 bp of any LRS call are deleted**, whatever the type or size (SV_core_LRS.cpp:2615-2626).
- **Surviving SRS calls are forced to 0/0 if >3 LRS reads span the site** (SV_core_LRS.cpp:2675-2676). With high-coverage LRS, SRS contributes almost nothing beyond LRS gaps.
- **Without `-k`, only SRS polishing is skipped** (ins_del.cpp:3192). A warning is printed (sv_call_main.hpp:259-261).
- **Even with `-k`, the polished contig is not written back.** `contig_polsihing_handler.run(bin_contig, ...)` edits a member buffer (ins_del.cpp:3206-3210). `LRS_SV_generating_Germline` regenerates `bin_contig` from the unpolished `target_path_idx_l[...].contig` (ins_del.cpp:3357, 3366). Polishing appears to have no effect on output. This is inferred from these lines; a test run comparing with and without `-k` would confirm.
- **9-mer profile format**: a header line, then `KMER<TAB>prob`. k is inferred from the first k-mer (RefHandler.hpp:340-361).
  - prob <0.6 or homopolymer marks the k-mer for polishing (RefHandler.hpp:386).
  - The loader allocates 4^k floats (RefHandler.hpp:354-355).
  - No generator ships in `src`, so a profile for your own ONT/HiFi data must be produced externally. Generator unverified.
- **ERR_PRONE preset** disables small-variant output (sv_call_main.hpp:218-220). It also changes SRS re-genotyping rules (SV_core.hpp:2535, 2542-2543).

### ngs_merging / ngs_fc

- **Neither tool reads usable CLI input.**
  - `ngs_merging` ignores argv and opens hardcoded `/media/fenghe/...` GRCh38 and 1KGP paths (NGS_merging.hpp:182-185, 223).
  - `ngs_fc` reads sites from hardcoded `/home/fenghe/E` (forceCalling.hpp:873).
- **`ngs_merging` uses an uninitialised bucket size.** `SV_MERGE_IDX_B_SIZE` is set only in `init()`, which `SV_merging` never calls (NGS_merging.hpp:163-166, 182-189).
- **`ngs_fc` skips all sites with pos < 23,504,400** (forceCalling.hpp:935). This is a debug leftover.
- **`ngs_fc` writes no results.** The GT is computed and discarded (forceCalling.hpp:979). Long INS/DEL and complex sites are never genotyped (983-997).

## Pitfalls D: population CSV, FC_BED, MEI-TR, tools

- **Documented FC_joint_ana command fails.** It uses lowercase `fc_joint_ana` and a 4th argument `1` (cohort_csv_analysis.md:241). The registered name is `FC_joint_ana` (src/SVcalling_core/ANALYSIS/analysis.cpp:167), and the 4th argument must be `HPRC` or `ONT`, otherwise xassert fires (FC_joint_calling.hpp:3003-3014). Tested: lowercase gives "wrong command".
- **Step 4 of the cohort doc is missing `-B PURE_STR`** (cohort_csv_analysis.md:167). `contig_file_split` parses `REF_REGION=` and `CONTIG_STR=` (analysis.hpp:3022, 3121), which only PURE_STR writes (SV_core_LRS.cpp:2398-2406). The default is VCF (sv_call_main.hpp:132).
- **Step 7.2 awk uses `$6,$7,$8` of a 3-column `bedtools merge` output** (cohort_csv_analysis.md:145, 248). Build prefixes from the `CONTIG_IN_REGION/` file names instead. Those names are BED ±500 bp (analysis.hpp:3026-3029; SV_core.hpp:2230).
- **`csv_region_500w` ignores every deletion.** The filter `SV_length < min_SV_LEN` (analysis.hpp:2999) runs on signed SVLEN, and HitSV writes DELs as negative (NovaSVRst.hpp:674; demo/lrs.vcf `SVLEN=-840`). Tested with two adjacent 60-bp DELs: no output.
- **`csv_region_500w` needs phased GTs.** It tests raw BCF bytes 4/5 (analysis.hpp:3008-3011; src/htslib/htslib/vcf.h:608-609). Unphased `0/1` is lost and `1/1` counts once. Keep `-f 1` (sv_call_main.hpp:133).
- **Hard cap of 24 sequences** (human chr1-22, X, Y):
  - csv_region_500w processes indices 0-23 and writes out of bounds for higher indices (analysis.hpp:2857, 2974, 3009).
  - FC_BED leaves regions on sequence index ≥24 uninitialised (SV_core.hpp:876-882).
  - entropy/small_var_den have the same problem (analysis.hpp:63-69, 536-542).
  - genomic_contest_* also cap at 24 (Genome_context_analysis.hpp:220, 251, 986).
  - For lizard/mole assemblies, put the target scaffolds first in the FASTA and BAM header, or rename and reorder them.
- **HOM detection depends on `_chr` in contig names** (FC_joint_calling.hpp:786-792, 813-826). With references named `1`, `NC_…` or `scaffold_…`, no contig is ever marked HOM. Haplotype counts in the sub-region stats (FC_joint_calling.hpp:2464) then undercount homozygotes. With `chr` names, a sample with one contig is called HOM, which is also wrong if the second haplotype simply failed to assemble.
- **FC_BED intervals that cross a 5-Mb segment boundary are silently dropped.** A region must be fully inside the current segment (SV_core.hpp:889; RefHandler.hpp:30, 246-248). For example, 4 999 000-5 001 000 is never assembled.
- **The "≤50 kb" limit is advisory.** No BED length check exists (SV_core.hpp:862-884). Contigs >50 kb are dropped in PURE_STR, but up to 3 Mb are kept with the ASM preset (SV_core_LRS.cpp:2324-2328, 2367). FC_joint_ana silently drops contigs >30 kb (FC_joint_calling.hpp:702-710), so 30 kb is the effective limit for population clustering.
- **Invalid `-B` gives empty output.** The README example `-B BED … -B PURE_STR` (README.md:365) uses a value that is not recognised (sv_call_main.hpp:237-242). With OUT_MODE_UNSET, the FC branch writes neither VCF nor PURE_STR (SV_core.hpp:2476-2479).
- **Read support is lost in pooled FASTA.** FC_joint_ana sets FULL_COVER=99 for every contig (FC_joint_calling.hpp:805), so `is_High_depth(0)` always passes (2783). Filter on SUP_R before `contig_file_split` (cohort step 5).
- **FASTA must be strict 2-line** (FC_joint_calling.hpp:798-805), and record 0 must be the ref (3033). Line-wrapped FASTA silently breaks pairing.
- **The FC_joint_ana MEI definition is substring-based**: `SINE`, `LINE`, `Retroposon`, `DNA`, `RC`, `LTR` (FC_joint_calling.hpp:172-180). The tandem-repeat definition is `Simple_repeat` or `Satellite` (165-169). A custom satellite library must use a class containing `Satellite` to be treated as VNTR. A SINE class must not contain `Satellite`.
- **The ONT vs HPRC flag changes only the "novel insert" threshold**: >60% vs >20% non-VNTR bases (FC_joint_calling.hpp:1316-1323). This decides whether a TR-rich SINE insertion counts as TR (same cluster) or novel (new cluster).
- **Minor bug in the nest flag.** The combined `is_many_mismatch` ORs `many_mismatch_2_ref` twice and omits `many_mismatch_3_ref` (FC_joint_calling.hpp:1777).
- **The README MEI-TR command is wrong.** It shows a 1-argument call (README.md:395), but the script needs `-f RM.out TRF.out` (mei_tr_array_single_region_detection.py:400-412). `-l` mode always raises TypeError (391).
- **The MEI-TR script crashes when a sample has no TRF block.** `calculate_coverage` returns `[]` (mei_tr_array_single_region_detection.py:78-80), and then the `:.2f` format raises (376). This happens before `mei_data.txt` is written (381). Tested.
- **MEI-TR pools all families of one class.** Grouping is by the first matching keyword, e.g. all `SINE/*` together (mei_tr_array_single_region_detection.py:211-215, 243). Strand and RM ID are ignored (197-200), so fragmented hits of one SINE count as several array units.
- **MEI-TR leftover variable.** `avg_distance_for_CV` (mei_tr_array_single_region_detection.py:306) is computed and never used. The spacing CV uses the raw mean (310), so arrays with near-zero spacing rely on the `stdev ≤ 15` branch (312).
- **The MEI-TR spacing list is padded with its own mean** (mei_tr_array_single_region_detection.py:263-265). DBSCAN distances in the chosen cluster are taken from the original adjacency and not recomputed between cluster members (292, 300).
- **The preprocessing script needs external header lines.** It expects `Sample <x>: File exists`, which HitSV never writes (contig_clustering_data_preprocessing.py:84). Its region regex is `chr[0-9XYM]+` (110). The default output filename contains Chinese characters (181).
- **The regulatory-annotation script is human-only.** It hardcodes GRCh38 Ensembl/GENCODE files in CWD (mei_tr_array_gene_regulatory_annotation.py:238-239). It splits region names on `_` (17-20), which breaks on `NC_…` or `scaffold_…` names.
- **Help strings are wrong for several tools.** `entropy_analysis`/`small_var_den` (analysis.cpp:173-174) differ from the real usage (analysis.hpp:455, 682). The region loaders for these tools skip the first BED line (analysis.hpp:61, 534). `align_seq` segfaults without its 9th argument (FC_joint_calling.hpp:2927, 2942; tested).
- **Developer paths left in code.** `/home/fenghe/E` is in `forceCalling_run`, reachable as `HitSV ngs_fc` (forceCalling.hpp:873; main.cpp:113). `/media/fenghe/...` is in unregistered `vcf_add_alt_string` (analysis.hpp:2243-2244).
- **FC_BED disables normal calling.** When `-b` is set, the regular germline pipeline (including SRS) is skipped (SV_core.hpp:2464, 2492). FC is ignored for tumour/normal runs (2464).
- **Demo-only issue.** With the 1-Mb demo reference, FC_BED aborted at RefHandler.hpp:265 ("region for building index must be within the current reference region") unless `-S 0 -E 0 -s 0 -F 1000000` was given (tested).
