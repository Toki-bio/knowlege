# Core data structures

Scope: HitSV v2.0.4, commit 93880367ef76f6440c918f7dd9744dda2659bcdb. All paths are relative to `src/`. Each line number was checked against the source.

## Ownership overview (who creates and who frees)

`NOVA_SV_handler::run` creates every top-level object (`SVcalling_core/sv_call_main.hpp:252-312`).

| Object | Created | Freed |
|---|---|---|
| `RefHandler refHandler` (stack) | `sv_call_main.hpp:253-255` | `destory()` (`RefHandler.hpp:215-218`) exists but `run` never calls it |
| `RefLocalRepeat *rlr` | `xcalloc`, only when `-n` is given (`sv_call_main.hpp:263-266`) | never freed |
| `BAM_STATUS *bs` | `xcalloc` (`sv_call_main.hpp:272`); filled from the `-T` JSON (`:275`) or computed on the fly (`:277`) | never freed |
| `KMER_ERROR_PROFILE_HANDLER` (stack) | `sv_call_main.hpp:279-282` | the `xcalloc`'d array (`RefHandler.hpp:355`) is never freed |
| `SV_CALLING_Handler *sve_h` | `new` (`sv_call_main.hpp:285`), then `init` (`:291-296`) | `distory()` (`:311`, defined at `SV_core.hpp:2335-2353`); no `delete` |

Main loop: `while(refHandler.load_seg_index())` builds the k-mer index (VCF mode only) and runs `SV_calling_in_region_main_pipeline()` (`sv_call_main.hpp:298-309`).

---

## 1. `RefHandler` (`SVcalling_core/RefHandler.hpp:157-329`)

Holds the reference and splits it into blocks.

| Field | Type | Meaning | Cite |
|---|---|---|---|
| `st_chr_ID, st_pos, ed_chr_ID, ed_pos` | int | region from `-S/-s/-E/-F` | 161-162, 196-199 |
| `referenceFilename` | char* | `-r` FASTA | 163, 200 |
| `c_ref_idx` | `faidx_t*` | FASTA index; `.fai` is built if missing | 165, 206; `clib/bam_file.c:262-274` |
| `curRefRegion` | RefRegion | chromosome currently held in `ref` (whole chromosome) | 168, 324 |
| `cur_r` | RefRegion | current block (always inside `curRefRegion`) | 169, 265 |
| `ref` | `{n,m,uint8_t*a}` | 2-bit/DNA5 copy of the **whole current chromosome** | 172, 314-321 |
| `ref_kmer_index` | `unordered_map<int32_t,RefIndexItem>` | 16-mer index of the block | 174, 273-291 |
| `bam_header_items` | vector<{name,len}> | copied from the BAM header | 176-186, 226-233 |
| `REF_IDX_KMER_LEN` | int | 16 | 209 |

Behaviour:
- `load_seg_index()` (243-267) moves to the next block of `SEGMENT_LEN` = 5,000,000 bp (`:30`). Blocks are contiguous and do not overlap (see io-formats.md §2). When the chromosome changes, `load_reference(chr_ID)` fetches the whole chromosome **by its BAM header name** (`:314`).
- `build_ref_index()` (273-291) indexes every 3rd 16-mer position (`i % 3 == 0`, `:283`) across the block. `search_ref_idx()` returns a position only for k-mers seen exactly once (`count == 1`, `:297`).
- Unused: `REF_BLOCK_LEN`, `REF_BLOCK_NUM`, `HASH_LEN`, `KMER_LEN` (`:26-29`) and `REF_index` / `REF_HASH` (`:33-68`). A grep finds no use of them outside this header.

## 2. `RefLocalRepeat` (`RefHandler.hpp:88-142`) and `RefLocalRepeatItem(R)` (`:70-86`)

Holds the `-I` "local repetitiveness" table.

| Field | Meaning | Cite |
|---|---|---|
| `repeatInfoMap` | `unordered_map<uint64_t, {int repeat_level}>`; key = `pos + (chrID << 32)` | 89, 120 |
| `region_length` | window length (600), read from the header line | 90, 113 |
| `region_step_size` | step (300), read from the header line; the in-class default of 300 is bypassed because the struct is `xcalloc`'d (`sv_call_main.hpp:264`) | 91, 114 |

- `load()` keeps only rows with `repeat_level >= 40` (`:123`).
- `search(chr, bg, ed)` scans every integer position in `[bg - bg%step, ed - ed%step + step]` and returns the hits (`:130-141`).
- `get_middle()` returns `position + len/2` and `range = len/2` (`:93-96`). This is how a window becomes a VNTR signal (`SV_core_SRS_ins_del.cpp:663-671`).
- `RefLocalRepeatItemR {chrID, position, repeat_level}` is the query-result record (`:74-86`).

## 3. `BAM_STATUS` (`SVcalling_core/bam_stat.hpp:14-218`) and `SIGNAL_PARAMETER` (`SV_core.hpp:44-83`)

| `BAM_STATUS` field | Meaning | Cite |
|---|---|---|
| `minInsertLen / middleInsertLen / maxInsertLen` | 1% / 50% / 99% quantiles of the insert size | 18-20, 126-128 |
| `isize_distribution` | fraction of read pairs for each isize in `[min, max)` | 21, 133-137 |
| `analysis_read_length` | modal read length (the first length bin holding >60% of reads) | 22, 105-118 |
| `ave_read_len`, `ave_read_depth` | mean read length; mean depth over sampled slices | 23-24, 116, 125 |

`SIGNAL_PARAMETER` is derived from it:
- `insert_size_min/max` and `MaxReadLen = analysis_read_length` (`SV_core.hpp:61-64`).
- `insert_region_len = max - 2*readlen` (`:65`).
- `max_del_dup_length = 50000` (`:67`).
- `reset_depth()` sets the depth-scaled cut-offs (`:70-81`), from the per-block depth: `SVE_MIN_SOLID_SCORE = max(2, 0.08*depth)`, `SVE_MIN_READ_NUM = max(3, 0.1*depth)` and `SVE_combine_min_score_step1 = clamp(0.4*depth, 16, 29)`.

## 4. `BAM_handler` / `READ_record` / `Ass_Block` (`SVcalling_core/ReadHandler.hpp`)

Stores the reads for one block.

| `BAM_handler` member | Meaning | Cite |
|---|---|---|
| `file`, `TL_file`, `re_alignment_file` | main BAM (`-n`/`-l`/`-X`/`-x`), the `-L` translocation BAM and the `-R` realigned BAM | 927-929, 613-627 |
| `read_list`, `um_read_list` | `vector<READ_record>` for mapped/SR reads and for unmapped reads | 921-922 |
| `base`, `quality`, `cigar` | concatenated kvec buffers; a `READ_record` keeps offsets into them | 946-948 |
| `storeReadBuff` | **300 MB** scratch buffer per handler (`new uint8_t[300000000]`) | 609 |
| `sr_read_index`, `um_read_index` | `SEGMENT_LEN/100 + 1` entries each | 606-608, 937-939 |
| `tr_loader` | `TransReadLoader` for discordant / translocation reads | 924, 623-626 |

- `READ_record` (24-41) fields: `flag`, `direction`, `position`, `read_l`, `read_end`, `cigar_l`, `read_index` / `quality_index` / `cigar_index` (offsets), `soft_left/right`, `NM_NUM`, `r_name`, `primary_chrID/pos`, `MAPQ`.
- `Ass_Block` (191-) is the input to local assembly: `ass_read_list` (per-read `ASS_reads_info` 88-109), `reads` (ACGTN strings) and the output `contigs`. Bases with quality below `MIN_BASE_QUAL` = 6 are masked to N (`:204-221`).
- `SV_CALLING_Handler` owns four instances: `SRS_read`, `LRS_read`, `LRS_Tumor_read`, `LRS_Normal_read` (`SV_core.hpp:2001-2005`). Each is opened in `init` (`:2262-2293`) and destroyed in `distory` (`:2336-2345`).

## 5. `SVE` / `SVE_core` / `SAMPLE_INFO` (`SVcalling_core/sve.hpp`)

The short-read (SRS) signal-cluster object.

| Item | Meaning | Cite |
|---|---|---|
| `SVE_core::r1, r2` | the two breakpoint intervals (RefRegion) | sve.hpp:36-37 |
| `SAMPLE_INFO::type_ID` | `SV::{DEL,DUP,INS,INV_1,INV_2,TRA,TRA_INV}` | 23, 147 |
| `signal_type` | `SIG::{DR,SH,...,LRS}`; LRS is printed as "TGS" | 19-20, 148 |
| `is_solid` | `RST::{BEGIN,END,SOLID,UNKNOWN}`: which side(s) have evidence | 28, 149 |
| `scoreB/E`, `readNumB/E` | `uint16_t` scores and read counts for begin/end | 153-156 |
| `vcfFilter()` | SOLID needs B+E ≥ 50 and each side ≥ 15; otherwise B+E ≥ 100 | 159-161, 219-227 |

- Storage: `SVE_L SRS_sve[SIG::LEN][SV::LEN]` plus `SRS_sveVNTR` (`SV_core.hpp:2038-2039`).
- The copy constructor is a raw `memcpy` (`sve.hpp:196-198`).
- `sve_combine_STEP1_SH()` merges signals whose r1 and r2 overlap and drops results with `score < min_score` and `< 2` reads (`sve.cpp:13-32`).

## 6. `LRS_SIG` / `LRS_SIG_Handler` (`SV_core.hpp:421-850`)

The long-read (LRS) signal object.

- `LRS_SIG {type, chrID, POS, END, isUsed, support_read_Num}` (`:421-428`).
- `LRS_SIG_Handler` holds `cigar_sig`, `clip_sig`, `small_sig`, `FINAL_sig` and an `LRS_SA_Handler` (`:830-850`).
- `LRS_SIG_COLLECTION` (`SV_core_LRS.cpp:1945-2178`) creates the signals:
  - CIGAR I/D of at least `MIN_sv_len` (`:2067-2076`, `:2083-2093`);
  - soft clips over 1000 bp (`:2042-2063`);
  - smaller indels as SMALL signals (ONT_Q20/HIFI only, `:1954-1955`, `:2077-2079`).
- The handler is a stack local of each pipeline call (`SV_core.hpp:2470`, `:2500`).

## 7. `NOVA_SV_FINAL_RST_item` + `NOVA_SV_TGS_INFO` + `TUMOR_FLAG` (`SVcalling_core/NovaSVRst.hpp`)

The final call; one item becomes one VCF record.

| Field | Meaning | Cite |
|---|---|---|
| `chr_ID, st_pos, SV_length, endPos` | position. `SV_length = alt.size() - ref.size()` (signed; INV uses the suggested length) | 1764-1766, 1874, 425-431 |
| `SV_type`, `bnd_type` | `"INS"`/`"DEL"`/`"INV"`/`"BND"`/`"INDEL"`/`"SNP"`; the BND subtype | 1871, 1875 |
| `ref, alt` | allele strings | 1864-1865 |
| `contig`, `cigar`, `contig_st_pos_in_region`, `region_ref_global_position` | assembled contig and its alignment | 1888-1893 |
| `signal_region_st/ed` | → `INFO/SIGNAL_REGION` | 1884-1885, 687 |
| `region_support_number[3][3]` | SRS read counts; `[2][0]` = ALT, `[2][1]` = REF, `[2][2]` = unknown | 1898, 610-612 |
| `final_genotype` | 0=0/0, 1=0/1, 2=1/1, 3=0\|0, 4=0\|1, 5=1\|0, 6=1\|1 | 1900-1903 |
| `LRS_INFO` (`NOVA_SV_TGS_INFO`) | `supp/not_supp/unknown_read_n`, `global_region_ID`, `haplotype_ID`, `total_hap_number_in_region`, `SV_in_HAP_ID`, `total_SV_number_in_HAP` | 290-328, 1862 |
| flags | `is_LRS_SV`, `is_BND`, `BND_BP_is_polishing`, `SVLEN_is_Imprecise`, `force_ALU`, `isVNTRcallerRst`, `QUAL_presetVNTR`, `passVNTR_depthFilter`, `deletion_pass_DR_filter`, `is_duplication_sv` | 1759-1762, 1904-1914 |
| `tumor_flag` | somatic class: TUMOR_UNIQ, NORMAL_UNIQ, TUMOR_NORMAL_SHARE or UNKNOWN | 330-359, 1863 |

- Result vectors in the handler: `result_SRS_INDEL`, `result_NGS_BND`, `result_BGS_INV` and `result_LRS_final` (`SV_core.hpp:2065-2069`).
- Created by `add_to_vector()` (`NovaSVRst.hpp:500-509`) and `store_BND()` (`:557-563`).
- Written by `writeVCF_final()` (`:622-746`).

## 8. `SV_CALLING_Handler` (`SVcalling_core/SV_core.hpp:1795-2653`)

The per-run engine; it owns every buffer.

| Group | Members | Cite |
|---|---|---|
| reference | `ref_handler`, `rlr`, `rlr_result`, `kmer_profile_handler` (non-owning pointers) | 1986-1992, 2250-2251 |
| BAM | file names + four `BAM_handler`s, `cur_header` | 1994-2007 |
| SRS signal data | `rdc` (depth counter), `DR/SH/UM_*_distribution`, `SRS_sve`, `SRS_sveVNTR`, `sig_para` | 2012, 2027-2039 |
| assembly | `ass_block`, `am` (`new MainAssemblyHandler[1]`), `dbgHandler` (`new MainDBGHandler1[1]`), aligners `ca`, `somatic_same_ca`, `ca_bnd`, `ca_re_locate`, `ga` | 2041-2056, 2072-2081, 2301-2307 |
| output | `vcf_w` (stdout or file), `vcfBuffStr` (1 MB), result vectors, `global_SV_ID` | 2061-2069, 2083, 2310-2315 |
| options | `MIN_sv_len`, `force_calling_ALU`, `random_phasing`, `output_small_var`, `preset_platform`, `Contig_aln_preset`, `OUT_MODE_int` | 2020-2026 |

Other details:
- `init()` (2237-2333) does the following:
  - opens the BAMs;
  - for SRS, builds the breakpoint distributions from `BAM_STATUS` (`:2269`);
  - allocates `alu_handler` (`:2275`);
  - registers the BAM header into `RefHandler` (`:2298`);
  - writes the VCF header unless `-H` or PURE_STR (`:2314-2315`).
- `distory()` frees `am` but not `dbgHandler`. It also calls `fclose(vcf_w)` even when `vcf_w == stdout` (`:2347`, `:2351`).
- `SV_calling_in_region_main_pipeline()` (2355-2637) has three exclusive branches:
  1. somatic, when `-X` is given (2360-2460);
  2. LRS forced-contig regions, when `-b` is given (2464-2490);
  3. basic germline (2492-2636). This runs the LRS calls (2495-2565), then the SRS calls (2567-2635). When both are present, the LRS results are re-genotyped with SRS reads (2523-2549) and the SRS calls are compared against LRS calls by `Hybrid_regenotype_NGS_USING_TGS_data` (2612-2615).

## 9. `KMER_ERROR_PROFILE_HANDLER` (`RefHandler.hpp:330-403`)

- The `-k` file is tab-separated `kmer\tprob` with one header line, which is skipped (`:348`).
- The k-mer length comes from the first data row, and the table holds 4^k floats (`:352-356`).
- `search_profile_and_indexing_region` marks k-mers for polishing (`'P'`) when they are all-A, all-T or have `prob < 0.6` (`:386`). It marks a 25-k-mer margin around them as edge (`'E'`) (`:364`, `:392-401`).
- `demo/9mer_profile` has the header `kmer\tprob` and 262,145 lines, i.e. 4^9 rows + header.

## 10. `RefRegion` (`cpp_lib/RefRegion.hpp`)

- Fields: `int chr_ID, st_pos, ed_pos` (`:135-137`), inclusive on both ends; `getLen() = ed - st + 1` (`:120`).
- The constructor and `set()` take `uint16_t chr_ID_` (`:25`, `:30`, `:39`). A contig index above 65535 wraps around.
- Both assert that the values are non-negative (`:27`, `:32`). `set_no_check()` skips the check and is used for the initial negative start (`RefHandler.hpp:202`).
