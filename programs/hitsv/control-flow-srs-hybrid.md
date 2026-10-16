---
title: "hitsv: Control flow: SRS and hybrid modes"
type: reference
generated: true
---

> Part of the [hitsv](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/hitsv/references/control-flow-srs-hybrid.md)

{% raw %}
# Control flow: SRS and hybrid modes

HitSV v2.0.4, commit 93880367ef76f6440c918f7dd9744dda2659bcdb. Paths are relative to `src/SVcalling_core/` unless prefixed. Every line number below was read directly from source.

## 0. Entry and per-block dispatch

1. `NOVA_SV_handler::run` (sv_call_main.hpp:252). With `-n` (SRS BAM): allocates `RefLocalRepeat` and loads `-I` (sv_call_main.hpp:263-266). `-I` is mandatory for SRS (check at sv_call_main.hpp:197-200), so the "empty genome context" warning branch (sv_call_main.hpp:267-270) is unreachable from the CLI.
2. BAM stats come from the `-T` JSON (sv_call_main.hpp:273-275). Without `-T` they are computed on the fly (sv_call_main.hpp:277; bam_stat.hpp:140-145).
3. The `-k` profile is loaded only when given (sv_call_main.hpp:280-282). If LRS and SRS are both given without `-k`, HitSV prints a warning only (sv_call_main.hpp:259-261).
4. `SV_CALLING_Handler::init` (SV_core.hpp:2237). For SRS it does the following:
   - `sig_para.init(bs)`: insert_size_min/max are the 1st and 99th percentiles of insert size (bam_stat.hpp:126-128), and `insert_region_len = max_isize - 2*readlen` (SV_core.hpp:60-68).
   - Builds the depth counter from `bs->ave_read_depth` (SV_core.hpp:2267).
   - Builds the breakpoint distributions (SV_core.hpp:2269; bam_stat.hpp:159-217).
   - Always constructs the ALU handler (SV_core.hpp:2275-2277).
5. The genome is processed in 5 Mb blocks (`SEGMENT_LEN 5000000`, RefHandler.hpp:30; loop at sv_call_main.hpp:298-309). Blocks follow the BAM header order (RefHandler.hpp:243-258).
6. Per block, germline and non-FC (SV_core.hpp:2492). If an LRS BAM is present, the LRS pipeline runs first (SV_core.hpp:2495-2565). The SRS pipeline follows (SV_core.hpp:2567-2635) in this order:
   1. `SRS_GET_SIGNALS_FROM_READs` (SV_core.hpp:2588)
   2. `reset_depth` (2591)
   3. `SRS_CLUSTERING_AND_COMBINE_signals` (2593)
   4. `SRS_SMALL_VAR_CALLING_PROCESS` (2599)
   5. `SRS_VNTR_VAR_CALLING_PROCESS` (2600)
   6. `SRS_LARGE_DELETION_CALLING_PROCESS` (2601)
   7. dedup and genotyping (2603-2604)
   8. INV (2606)
   9. BND (2608)
   10. hybrid merge, only if LRS is present (2612-2615)
   11. write (2624-2634)
   - Note: loading realigned reads from `-R` is disabled with `if(false)` (SV_core.hpp:2590). The `-R` option has no effect on calling.

## 1. SRS signal types, detection and thresholds

### 1a. Read-level filters (SV_core_SRS_ins_del.cpp:419-474)

| Step | Rule | Lines |
|---|---|---|
| Duplicates skipped; depth counter incremented | before the secondary/supplementary filter, so those reads count toward depth | ins_del.cpp:437-441 |
| Secondary/supplementary skipped | | ins_del.cpp:444-447 |
| Unmapped reads go to the UM list | after the compact filter | ins_del.cpp:449-452, 365-375 |
| MAPQ filter | MAPQ>=20, OR mate MAPQ>=20 and same chr and \|isize\|<20000 | ins_del.cpp:17-25, 453 |
| Long-range/inter-chr pairs kept for TL loading | `mtid!=tid` or \|pos-mpos\|>200000 | ins_del.cpp:471-473 |
| Block depth | `read_counter*MaxReadLen/blockLen` | ins_del.cpp:476 |

### 1b. Signal detectors

| Paper term | Code | Detection rule | Signal score / region | Lines |
|---|---|---|---|---|
| CLP (clipping) | `SRS_handleSASignal` -> `SRS_storeClipSignals`, stored as `SRS_sve[SH][INS]` | Read must have a clip or `max(gaps,NM)>=2` (`bam_aligned_analysis`). It must pass `storeReadSR` filters: compact-length, a poly-A/T clip-tail filter, and a low-quality clip filter. Any soft clip >=1 bp (`MIN_CLIP_LEN 1`) emits a signal. | score `MIN(15,MAPQ)`; region = clip point +/-10 bp | ins_del.cpp:393-417, 457-463, 176-188, 51-62; ReadHandler.cpp:78-92 |
| SDA (CIGAR divergence) | `SRS_storeMismatchSignals`, stored as `SRS_sve[SH][INS]` | Requires NM>=3 and no clipping (ins_del.cpp:189-191). Counts mismatches with base quality >=20, plus every I/D base (ins_del.cpp:64, 86-98). Requires >=3 independent (non-adjacent) events (ins_del.cpp:107-114). | score `MIN(3,MAPQ)`; region = 200-bp bin `pos - pos%200` .. +200 | ins_del.cpp:157, 166, 173 |
| DRP (discordant pairs) | `bam_is_DR_signal` -> `SRS_handleDRSignal`, stored as `SRS_sve[DR][type]` | Pair is discordant if any of: mates on different chromosomes; same orientation; isize<=1%ile or >=99%ile; wrong FR order. Type assignment: DEL (isize > max), DUP (isize < min or RF orientation), INV_1/INV_2 (same orientation), TRA/TRA_INV (different chromosome). Pairs with \|isize\| >= `max_del_dup_length` (50000) fall to no type. | score `MIN(15,MAPQ)`; region length = `insert_region_len` | clib/bam_file.c:509-530; ins_del.cpp:540-563; sve.hpp:168-185; SV_core.hpp:67 |

Notes:
- DUP signals are generated (ins_del.cpp:550,552) and clustered (ins_del.cpp:626-639), but no caller consumes `SRS_sve[DR][DUP]`. A grep for `SV::DUP` finds only these two lines.
- `SRS_sve[DR][INS]` is printed (SV_core_SRS.cpp:197) but never assigned by `SRS_handleDRSignal`.

### 1c. Clustering and thresholds (`SRS_CLUSTERING_AND_COMBINE_signals`, ins_del.cpp:622-677)

1. Step 1, `SRS_single_type_sve_combine` (ins_del.cpp:239-362):
   - Overlapping signals of the same begin/end side are pooled.
   - A cluster needs **>=3 signals** (`MIN_ACC_READ=3`, ins_del.cpp:251,257).
   - The summed score must exceed `SVE_combine_min_score_step1 = clamp(0.4*depth, 16, 29)` (SV_core.hpp:77-79; ins_del.cpp:315-317).
   - Both sides must span less than 5000 bp (ins_del.cpp:317).
   - The breakpoint is the argmax of the summed insert-size-derived (DR) or quadratic 10-bp (SH) breakpoint probability (ins_del.cpp:196-236; bam_stat.hpp:163-195). It must reach `2*max(single-read distribution)` (SV_core_SRS.cpp:349-358; ins_del.cpp:334-336).
2. Step 2, `SRS_sve_begin_end_combine` (ins_del.cpp:568-600): BEGIN and END sides are paired into SOLID.
   - A SOLID cluster must have score >= `SVE_MIN_SOLID_SCORE` and read count >= `SVE_MIN_READ_NUM`. A one-sided cluster needs **2x** both (ins_del.cpp:586-587).
   - `SVE_MIN_SOLID_SCORE = max(2, 0.08*depth)` and `SVE_MIN_READ_NUM = max(3, 0.10*depth)` (SV_core.hpp:72-75).
   - DUP/INV/TRA DR clusters use 1.5x these values (ins_del.cpp:633-636).
3. Repeat windows are added as VNTR candidates (score 99, always SOLID) (ins_del.cpp:659-671). See section 2.
4. Assembly gate `SRS_SV_region_depth_filter` (ins_del.cpp:901-916): skip if the max 32-bp-bin depth is >15x expected, or the mean is >5x expected. Depth is normalised by `ave_depth*32/readlen` (SV_core.hpp:155-159, 187-205). This gate is applied to every SRS candidate type (SV_core_SRS.cpp:268, 445, 485, 521, 577).

### 1d. Paper vs code

| Paper claim | Code | Verdict |
|---|---|---|
| 600-bp windows, 300-bp overlap | Only in `srs_fa_stat` repeat indexing (`REGION_LEN=600`, `REGION_STEP_LEN=300`, fa_stat.hpp:34-35). Signal discovery is overlap clustering inside 5 Mb blocks (RefHandler.hpp:30); there are no 600-bp windows. | **Partly**: true for the R index only |
| Candidate if R>40 | Kept if `repeat_level >= 40` (RefHandler.hpp:123). Each kept window becomes a candidate (ins_del.cpp:663-671). R takes values in {20,40,...,180} (fa_stat.hpp:36-38, 70). | **Confirmed**, with `>=` instead of `>` |
| SDA/CLP/DRP count > max(10% depth, 2) | Read count >= `max(3, 0.10*depth)` and score >= `max(2, 0.08*depth)` (SV_core.hpp:72-75). Doubled for one-sided clusters (ins_del.cpp:586). Also >=3 signals and score > clamp(0.4*depth,16,29) earlier (ins_del.cpp:315-317). | **Refuted in detail**: the floor is 3, not 2, and there are extra gates |
| Depth basis | Block-local depth from read count (ins_del.cpp:476, SV_core.hpp:2591), not the JSON depth | note |

## 2. Local repetitiveness R (`-I` file)

- **File format** (RefHandler.hpp:102-127): the header line is `region_len  step  0` (written at fa_stat.hpp:49). Each body line is `chrIdx  windowStart  R`, where chrIdx is the FASTA order index (fa_stat.hpp:52, 97). R is the largest k in {20,40,...,180} at which the 600-bp window still contains a repeated k-mer. `withloop` is never reset between k values (fa_stat.hpp:68-92), so the scan stops only when no repeat exists at the first k tried.
- **Load cutoff**: windows with R<40 are discarded at load time (RefHandler.hpp:123-125).
- **Decisions that depend on R:**
  1. Every R>=40 window overlapping the block becomes an `SRS_sveVNTR` candidate centred on the window (`get_middle`, RefHandler.hpp:93-96; ins_del.cpp:663-671). These candidates are assembled by the dBG VNTR caller (SV_core_SRS.cpp:563-590).
  2. In `SRS_SMALL_VAR_CALLING_PROCESS` the repeat lookup runs (SV_core_SRS.cpp:539-541), but the result is only printed (545-546). The decision to call VNTR uses `withRepeat` from assembly (SV_core_SRS.cpp:543), not R.
  3. Hybrid: when an LRS SV overlaps an R window and depth<10 and preset is `ERR_PRONE`, SRS re-genotyping of that LRS SV is skipped (SV_core.hpp:2528-2536).
- Search granularity: `search` scans every integer position from `bg - bg%step` to `ed - ed%step + step` and hash-looks-up each (RefHandler.hpp:130-141).

## 3. Insertion/deletion calling

### 3a. Path A: SH/INS signals -> greedy overlap assembler ("MEI_AND_LONG_SV")

1. Called from SV_core_SRS.cpp:507-561 with `call_ALU_TL=true` (SV_core_SRS.cpp:532). An unused TL-count gate `MIN_TL_read_number = max(4, 0.15*depth)` is computed (526) but `cur_TL_read_num` is 0 and the test is `>=0`, so it always passes (523, 529).
2. Read loading `SRS_assembly_load_read` (ins_del.cpp:1362-1409):
   - SR reads in main/supp regions, up to 1000 (1370, 1383, 1387)
   - UM reads +/-1000 bp (1391-1394)
   - TL reads (1402)
3. Assembly: `ass_block.run_assembly(am)` (ins_del.cpp:2025) uses the greedy/iterative word assembler. Word lengths are 26,31,...,121 (GreedyAssembler.hpp:217-219); minCoverage 2, minSupportReads 3, maxAssemblyCount 5 in repeat mode and 10 in normal mode (GreedyAssembler.hpp:221-225, 241-248). If the max word length is reached, `withRepeat=true` (GreedyAssembler.cpp:716-720; ins_del.cpp:2026).
4. Contig filters: drop contigs with a non-zero ending reason, length <100, or <=2 new support reads at k<100 (ins_del.cpp:2043-2048).
5. ksw2 realignment of the contig to the local reference: `SRS_alignment_and_get_var` (ins_del.cpp:866-899) -> `ca.align_non_splice_default_ref` (ksw_extd2_sse, Contig_aligner.hpp:160/176/206/267).
   - Contigs with NM>5% of matched length are rejected (ins_del.cpp:892-895).
   - Variants are extracted by `get_canditate_SVs` (ins_del.cpp:897).
6. If an SV was found, `SRS_read_cover_repeat_filter` may set `withRepeat` (ins_del.cpp:2103-2105). A true `withRepeat` triggers path B on the same signal (SV_core_SRS.cpp:543-553).
7. Large-insertion reconstruction (only when `call_ALU_TL`, ins_del.cpp:2110-2121):
   - **Poly-A/T tail joining** (`SRS_combine_repeat_tail_of_contigs`, combine_contigs.cpp:197-327):
     - Aborted if the local reference has >=2 homopolymer runs (>=13 bp) of the same base (combine_contigs.cpp:205-226).
     - Finds one contig with a homopolymer >=15 bp at its tail (within 80 bp of the end) and another with one at its head (within the first 80 bp), using the same base (combine_contigs.cpp:54, 76-86, 245).
     - Concatenates them through the longer homopolymer and realigns with ksw2 (combine_contigs.cpp:253-320).
   - **Transposition-based insertion** (`SRS_handle_tran_based_ins`, combine_contigs.cpp:569-721), run only if nothing has been called yet (ins_del.cpp:2114):
     - Contigs whose read-type runs are exactly SR+TL or TL+SR (combine_contigs.cpp:90-141) yield source breakpoints from TL reads.
     - With exactly one left and one right breakpoint (combine_contigs.cpp:605), the source reference segment between them is spliced in as the insert sequence (combine_contigs.cpp:623-653). Load windows: 5 kb for segments <5 kb, 50 kb for <50 kb (696-701).
     - The call is realigned and flagged IMPRECISE (combine_contigs.cpp:703-704, 716-720).
   - **Force ALU** (only if still nothing called and `force_calling_ALU`): see 3c.
8. Output: `resultFilter` then keep `will_be_output_to_vcf` (ins_del.cpp:2124-2136). INS/DEL are genotyped later (SV_core_SRS.cpp:217-234).

### 3b. Path B: VNTR / repeat -> coverage-guided dBG (`SRS_assembly_variations_INS_DEL_VNTR`, ins_del.cpp:2182-2334)

Path B runs for R-window candidates (SV_core_SRS.cpp:563-590) and for path-A signals with `withRepeat` (SV_core_SRS.cpp:550).

1. Reference: region +/-250 bp. Supp is joined when <10 kb away (ins_del.cpp:2187, 2209-2217).
2. Reads: +/-250 bp, minQUAL 3, up to 60000 reads (ins_del.cpp:2188, 2238).
3. k-mer error correction at k=65 (ins_del.cpp:2258-2262).
4. dBG at k = 60, 80, 100, 120 (ins_del.cpp:2270-2275), using k-mers with depth >=3 (`minDepthNGS=3`, ins_del.cpp:2269, 2283).
   - Higher k is tried only while complex nodes remain (ins_del.cpp:2304-2310).
5. **Traversal** `getDBGPathM1` -> `D_searchPath` (DBGAssembler.hpp:3321-3367, 3261-3319):
   - DFS from each key node to the stop set; a stop node with read support completes a path.
   - Per-node max visit count `MAX_nodePassTime = 13 - 4*complexMode`, i.e. **13, 9, 5, 1**, retried when too complex (DBGAssembler.hpp:3339-3347, 3288-3289).
   - Abort limits: >30000 node visits, >120 paths, or depth >100 (DBGAssembler.hpp:3264-3269).
   - **Paper "9" = the second tier only; the first attempt uses 13.**
   - Paths are then read-path filtered (`isPathFilter`, DBGAssembler.hpp:3378).
6. **Haplotype pair selection** (`analysisAndGenerateSV_M2`, DBGAssembler.hpp:1917-2038):
   - Every pair (i, j>=i), including homozygous i==j, is scored by `calFitRateM2` (DBGAssembler.hpp:1950-1963).
   - The score is the mean of `|preDepth*curCount - preCount*curDepth|/60` over unitig junctions, i.e. depth consistency (DBGAssembler.hpp:1543-1565 area; formula at line 1556).
   - The score is multiplied by 0.7 for equal lengths, 0.6 for 0/0, 0.85 for 0/1 (1635-1646), by 4 for >25 (or >15 when homozygous) junctions (1647-1648), then +0.03 per unused read-path (1649).
   - Lowest score wins (1964).
   - Reject if best > 0.60 / 0.40 / 0.30 when there are 1 / <=3 / more pairs (1983-1998).
   - QUAL 0-30 comes from the ratio to the first dissimilar alternative (2023-2037).
   - This is a heuristic depth-fit score, **not** a probabilistic likelihood.
7. Results are kept only if the set does not contain both INS and DEL (ins_del.cpp:2314-2332).

### 3c. MEI / ALU-specific logic (`force_calling_ALU`)

- `force_calling_ALU = true` is hard-set and has no CLI option (sv_call_main.hpp:175). It is passed to `init` (sv_call_main.hpp:294; SV_core.hpp:2277).
- `ALU_force_handler::init` (SV_core_SRS_inv_tra.cpp:820-859):
  - Embeds one **human Alu consensus**, 320 chars including N positions (inv_tra.cpp:822).
  - Defines 5 anchor k-mers per strand at positions 28/68/109/190/213 with lengths 20/18/16/17/16 (inv_tra.cpp:844-850).
- `SRS_force_calling_MEI_ins` (ins_del.cpp:1774-1971) runs on path A only, when no other SV was produced (ins_del.cpp:2119-2120):
  1. Exact substring search for an anchor in each contig (inv_tra.cpp:901-911).
  2. A synthetic insert of length `320 + (pos%15 - 7)`, padded with N (ins_del.cpp:1801, 1810-1838).
  3. A realignment check that requires an insertion of 280-350 bp in the CIGAR (ins_del.cpp:1866-1884; Contig_aligner.hpp:687-703).
  4. An Alu-consensus check on 200 bp: <=50 N and fewer than 10 unmatched bases (inv_tra.cpp:861-899).
  5. The position must be >=100 bp inside the signal region (ins_del.cpp:1940-1948).
  6. Emitted as INS, IMPRECISE, with the internal `force_ALU` flag set, which is not written to the VCF (ins_del.cpp:1960-1965).
  7. In genotyping, the REF-support count is halved for force_ALU calls (NovaSVRst.hpp:1577-1578).
- No other MEI family (L1, SVA, or non-human SINEs) is modelled. Generic MEIs depend on the poly-A/T joining and TL splicing paths (3a.7).

### 3d. Large deletions (DR/DEL)

`SRS_LARGE_DELETION_CALLING_PROCESS` (SV_core_SRS.cpp:473-505) calls path A with `call_ALU_TL=false` (501).
- The load extension grows with span: 2 kb by default, 5 kb when the span is >1800 and <5000, 50 kb when <50000 (489-495).
- Spans >=50000 are skipped (496-497).

## 4. Inversions and translocations

### INV (SV_core_SRS.cpp:592-605)

1. `SRS_INV_signal_combine`: DR INV_1 and INV_2 clusters with overlapping r1 and r2 are merged (236-265).
2. Depth gate, then a skip when span >50000 (268, 286-290).
3. `SRS_assembly_and_get_breakpoints_INV` (inv_tra.cpp:636-817):
   - Reads +/-readlen around both sides, up to 1000 (652-669).
   - Greedy assembly (674).
   - Contigs are aligned against 4 combined references built from +/-400 bp at each side in both orientations (677-700+).
   - Rejected if the two inversion junction lengths differ >1.5x or by >150 bp (793-801).
   - Breakpoints are the mean of the supported combinations (803-814).
4. GT from pair counts on both sides (SV_core_SRS.cpp:301-322):
   - 1/1 if alt > 4*ref; 0/0 if alt*4 < ref.
   - 0/0 if unknown > 3*(alt+ref) or alt < 2.
5. Stored if 30 <= len <= 50000. Explicit REF/revcomp ALT for <=2000 bp, otherwise `<INV_REF>/<INV_ALT>` (373-401).
6. Dropped if it overlaps an INS/DEL >=50 bp (330-341).

### TRA (SV_core_SRS.cpp:409-453)

1. Only `DR/TRA` and `DR/TRA_INV` are processed (419-439), with a depth gate (445).
2. `SRS_assembly_and_genetyping_BND` (inv_tra.cpp:11-210) **does no assembly**. Despite its name, it counts BND-supporting pairs vs reference-spanning pairs around the cluster midpoint (inv_tra.cpp:76-145).
   - GT: 2 if BND >= 3*REF; 0 if BND*6 <= REF; 0 if BND < 3 (inv_tra.cpp:152-158).
3. The record is stored with `BND_BP_is_polishing=true` (inv_tra.cpp:202). `is_BND_IMPRECISE()` is therefore always false (NovaSVRst.hpp:1198-1199), so every TRA with GT>0 is written (SV_core_SRS.cpp:458-463).
4. `SRS_assembly_and_get_breakpoints_TRA` exists (inv_tra.cpp:360) but has no call site; the only other match is its declaration (SV_core.hpp:2149).

### TL BAM (`-L`) use

- Opened in `BAM_handler::init` (ReadHandler.hpp:621-623).
- `TransReadLoader::load_read` uses the TL BAM if present (M1). Otherwise it scans the main BAM at each mate cluster (M2) (ReadHandler.hpp:388-389).
  - M1 matches TL records by (pos,mpos) keys from `trans_list` after swapping pos/mpos back (ReadHandler.hpp:405-433).
  - M2 clusters mates within 20 kb, skips low-complexity (>80% one base) or >3000-read regions, and requires >3 trans reads (ReadHandler.hpp:475-530).
- TL reads feed only the assembly read sets: path A (ins_del.cpp:1402) and the transposition-insertion reconstruction (combine_contigs.cpp:144-195, 569-721). **They are not used for TRA/BND calling.**

## 5. Hybrid (LRS + SRS)

Order per block (SV_core.hpp:2495-2635):

1. **LRS calls first** (2502-2520). They are then **re-genotyped with SRS reads** (2523-2549).
   - Skip when the LRS SV is in an R window and LRS depth <10 and preset is ERR_PRONE (2528-2536).
   - Skip SVs shorter than `MIN_sv_len` (2537-2538).
   - ERR_PRONE: keep the LRS GT unless REF is "perfectly" supported by SRS (2542-2543; NovaSVRst.hpp:615-620).
   - Set 0/0 if LRS has exactly 1 supporting read and >=1 non-supporting read, and the SRS non-supporting count is >=20x the supporting count (2545-2546; NovaSVRst.hpp:604-613).
   - Adopt the SRS GT if LRS has 0 non-supporting reads and <4 supporting reads (2547-2548; NovaSVRst.hpp:599-602).
2. **The full SRS pipeline runs independently** (2585-2609). The `-Y` INNER_TEST path instead loads SRS calls from a VCF (2579-2583).
3. `Hybrid_regenotype_NGS_USING_TGS_data` (SV_core_LRS.cpp:2636-2687):
   - Drops every SRS call within 300 bp of any LRS call on the same chromosome (SV_core_LRS.cpp:2607-2634).
   - Counts LRS reads fully spanning +/-500 bp of each remaining SRS call (2653-2673).
   - Sets the SRS GT to **0/0 if there are >3 spanning LRS reads**, and to 0/1 if there are 1-3 (2675-2677). VNTR SRS calls with only non-supporting LRS reads are set to 0/0 (2680-2682).
   - Net effect: **LRS calls are the backbone; SRS calls survive only where LRS coverage is <=3.**
4. **SRS polishing of LRS contigs** (ins_del.cpp:3190-3213): runs only when an SRS BAM is present **and** the `-k` profile is loaded (`kmer_error_profile != NULL`).
   - Each selected haplotype contig longer than 50 bp is polished against SRS reads in the region +/-300 bp (3199-3210).
   - `Contig_Polishing::run` (Contig_polishing.hpp:508-517) does three things:
     1. Indexes 13-mers in contig positions flagged by the profile (Contig_polishing.hpp:195, 202-224).
     2. Anchors SRS reads that have >25 unique-k-mer hits (243-260).
     3. Fixes 1-2-base differences where >=3 reads and >=70% agree (116-127, 476-495), and homopolymer lengths with >=3 reads and >70% agreement (378-410, 420-474).
   - **Caveat:** the polished `bin_contig` is not written back to `target_path_idx_l[...].contig`. `LRS_SV_generating_Germline` rebuilds `bin_contig` from the unpolished string (ins_del.cpp:3206-3217 vs 3357, 3366). Polishing therefore appears not to affect output (inferred from these lines; no write-back found).
   - Separately, ONT_Q20/CCS presets apply LRS-only homopolymer correction (ins_del.cpp:3170-3188).
5. **9-mer error-profile format** (`KMER_ERROR_PROFILE_HANDLER::load_data`, RefHandler.hpp:340-365):
   - TSV with a header line (`kmer\tprob`, skipped at RefHandler.hpp:348). Then `KMER<TAB>prob`.
   - k is taken from the first k-mer's length (RefHandler.hpp:352-356). The demo has a header plus 4^9 = 262144 rows (`demo/9mer_profile`, 262145 lines).
   - `prob` behaves like a per-k-mer accuracy. K-mers with prob <0.6 or all-A/all-T are flagged `P` (polish), with +/-25 k-mer `E` edges (RefHandler.hpp:364, 386-400).
   - No generator for this file was found in the source tree (grep of src for `9mer`/`KMER_PROFILE` hits only the loader/consumers). Its origin is unverified.
6. **Without `-k`**: the warning at sv_call_main.hpp:259-261 is printed, the profile stays NULL (RefHandler.hpp:333), and SRS polishing is skipped (ins_del.cpp:3192). Every other hybrid step (1-3) still runs.

## 7. Human-specific / non-human-relevant hardcoding (call path)

See pitfalls-C.md for the full list. Key items:
- Alu consensus and force-ALU always on (inv_tra.cpp:822; sv_call_main.hpp:175).
- `asm5` is the default contig-alignment preset, labelled "Recommended for Human"; `asm10`/`asm20` are recommended for other species (sv_call_main.hpp:150-153).
- The `-I` index uses FASTA order indices (fa_stat.hpp:52, 97), while lookups use BAM tids (RefHandler.hpp:130-136; ins_del.cpp:664). The FASTA and BAM header orders must therefore be identical.
- Insert-size and read-length values in `call` come from data (bam_stat.hpp:72-138); none are hardcoded there. `ngs_fc` hardcodes them (forceCalling.hpp:709-710, 889).
{% endraw %}
