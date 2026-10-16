---
title: "hitsv: Population CSV analysis, FC_BED contigs and MEI-TR arrays"
type: reference
generated: true
---

> Part of the [hitsv](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/hitsv/references/population-csv-mei-tr.md)

{% raw %}
# Population CSV analysis, FC_BED contigs and MEI-TR arrays

Source: HitSV at commit 93880367ef76f6440c918f7dd9744dda2659bcdb (v2.0.4). Paths are relative to the repo root. `src/htslib` is vendored and is cited only where the behaviour of an htslib call matters. Lines marked "tested" were also run against the locally built binary (built from source at the pinned commit, see `build-run-test.md`) or `python3` on synthetic or demo inputs.

## 1. End-to-end cohort CSV workflow

The README gives the outline in five steps (README.md:374-382). The full recipe is in cohort_csv_analysis.md. Each step below lists what the code actually requires.

| # | Step | Command (as documented) | Input | Output | Notes verified in code |
|---|------|------------------------|-------|--------|------------------------|
| 1 | Per-sample SV calling | `HitSVL call -p $PRESET -l $BAM -r GRCh38.fa -o D.vcf` (cohort_csv_analysis.md:81) | LRS BAM | `RST/<sample>/D.vcf` | Batch line in cohort_csv_analysis.md:86 leaves out the script name (`sbatch %s ...`). The binary is `HitSV`, not `HitSVL` (main.cpp:115; the built binary prints `Program: HitSV`). |
| 2a | Per-sample CSV windows | `HitSV tools csv_region_500w D.vcf 500 50 > S.main.bed` (cohort_csv_analysis.md:111) | the sample VCF | BED-like rows `chr st ed HAP1/HAP2 TRUE_POS/500W` (analysis.hpp:2873, 2878) | See section 1.1 for the algorithm and its limits. |
| 2b | Keep 500-bp windows per haplotype | `grep HAP1 \| grep 500W \| sort \| uniq` (cohort_csv_analysis.md:113-114) | S.main.bed | S.hap1.bed, S.hap2.bed | "500W" means 500-bp bins: `(pos/500)*500`, step 500 (analysis.hpp:2875-2879). |
| 2c | Count haplotypes per window across the cohort | `cat */*hap*bed ... \| sort \| uniq -c` (cohort_csv_analysis.md:129-130) | all hap BEDs | ALL_single_hap_summary.txt | The count is the number of haplotypes with a CSV in that window. |
| 3a | Frequency filter | awk `$1 >= 56` (cohort_csv_analysis.md:138) | summary | AF_0.05_CSV_AF.bed | The 56 is a cohort-specific allele count (labelled AF 0.05), not a code constant. |
| 3b | Sort and merge | `bedtools sort`, `bedtools merge` (cohort_csv_analysis.md:144-145) | AF BED | `*.sort.merge.bed` (3 columns) | Merged intervals must stay ≤50 kb (README.md:370). In practice FC_joint_ana skips contigs >30 kb (FC_joint_calling.hpp:702-710). |
| 4 | Force-generate local contigs per sample | `HitSVL call -p $PRESET -b merge.bed -l $BAM -r GRCh38.fa -o HitSV_FC_CSV_RST.txt` (cohort_csv_analysis.md:167) | BAM + BED | `RST/<sample>/HitSV_FC_CSV_RST.txt` | **You must add `-B PURE_STR`.** Steps 5 and 6 parse PURE_STR fields (`REF_REGION=` prefix skipped by `+11`, analysis.hpp:3022; `CONTIG_STR=`, analysis.hpp:3121). The documented command uses the VCF default (sv_call_main.hpp:132). See section 4 for the format. |
| 5 | Filter by read support | awk `split($4,A,"="); split(A[2],B,":"); B[1] > 1` (cohort_csv_analysis.md:192) | PURE_STR txt | `HitSV_FC_CSV_RST.FILTER.txt` | Field 4 is `SUP_R:...=n:...` (SV_core_LRS.cpp:2402). This keeps contigs with SUP_R ≥ 2 and drops STAT≠PASS lines. |
| 6 | Pool contigs per region | `HitSVL tools contig_file_split list.txt ./CONTIG_IN_REGION/ GRCh38.fa` (cohort_csv_analysis.md:207, 212) | list of FILTER.txt paths | One FASTA per region named `<chr>_<st>_<ed>`: first record `>ref<chr>_<st>_<ed>`, then `><sample>_<region>_<CI>` | The sample name is the **parent directory** of each listed file (analysis.hpp:3077-3078). Reference is fetched with faidx (analysis.hpp:3104-3109). CI comes from field 3 split on ":" at index 9 (analysis.hpp:3124-3126). Each region file is BED ±500 bp: no-contig lines get ±500 (analysis.hpp:3026-3029), and PASS lines carry the 500-bp load margin (SV_core.hpp:2230; tested: BED 869300-870500 gave `1_868800_871000`). |
| 7a | Annotate repeats | `RepeatMasker $P.fa`, `trf $P.fa 2 7 7 80 10 10 500 -h -d -ngs > $P.fa.trf.out` (cohort_csv_analysis.md:237-239) | region FASTA | `$P.fa.out`, `$P.fa.trf.out` | RepeatMasker is called with no `-species` or `-lib`. TRF must use `-ngs`, because the loader expects `@name` header lines (FC_joint_calling.hpp:853-856). |
| 7b | Cluster haplotypes | doc: `HitSVL call tools fc_joint_ana $P.fa $P.fa.out $P.fa.trf.out 1` (cohort_csv_analysis.md:226, 241) | the three files | 15 log files in **CWD** (see 1.3) | **As written, this command fails.** The registered name is `FC_joint_ana` (analysis.cpp:167), and lowercase is rejected (tested: "wrong command"). The 4th arg must be `HPRC` or `ONT`; any other value hits xassert (FC_joint_calling.hpp:3003-3014). Working form: `HitSV tools FC_joint_ana P.fa P.fa.out P.fa.trf.out ONT`. |
| 7c | Batch | awk on `$6,$7,$8` of merge.bed (cohort_csv_analysis.md:248) | merge.bed | — | `bedtools merge` writes only 3 columns, so `$6-$8` are empty. Build the prefixes from `ls CONTIG_IN_REGION/` instead. |
| 8 | Post-process | `python3 py/contig_clustering_data_preprocessing.py <file> [--output o]` (contig_clustering_data_preprocessing.py:168-185) | concatenated `sub_region_detail.log` content | `<input>.预处理1.txt` by default (line 181) | Expects `Sample <name>: File exists` header lines (line 84). HitSV never prints these (grep finds that string only in this .py file), so a user wrapper has to add them. |
| 9 | Optional region annotation | `python3 py/mei_tr_array_gene_regulatory_annotation.py regions.txt` (lines 225-236) | lines `chr_st_ed` (lines 16-21) | `_annotation_overlap.csv/.png`, `_annotation_details.tsv` (lines 210-212) | Hardcoded human GFF3/GTF file names (lines 238-239). |

**Where `nestedStructures.log` comes from:** `FC_joint_calling_handler::contig_clustering` writes it at FC_joint_calling.hpp:2911-2913, from `findNestedStructures` (FC_joint_calling.hpp:2699-2747), into `work_dir = "./"` (FC_joint_calling.hpp:2997). It is produced by `HitSV tools FC_joint_ana` (step 7b), in the directory you run it from (tested). A record is "nested" when `base_is_also_INS_in_ref` is true (FC_joint_calling.hpp:2711): the SV's base breakpoint lies inside an insertion of the base contig relative to the reference, i.e. `mapping_type == 1` (FC_joint_calling.hpp:1731). The record's `nest_string` has mismatch-density checks for contig3→2, 2→ref and 3→ref (FC_joint_calling.hpp:1758-1793).

### 1.1 `csv_region_500w` algorithm (analysis.hpp:2947-3018)

1. Reads only PASS records (analysis.hpp:2985-2986). Tested: LOW_DEPTH records gave no output.
2. Drops records with `SVLEN < min_SV_len` (analysis.hpp:2999). HitSV writes DELs with negative SVLEN (demo/lrs.vcf: `SVLEN=-840`; NovaSVRst.hpp:674), so **all deletions are dropped**. Tested: two adjacent DELs gave no output. INS end is set to POS+1 (analysis.hpp:3001-3002).
3. Haplotype assignment compares raw BCF GT bytes: `GT[0][0]==4` → HAP1, `GT[0][1]==5` → HAP2 (analysis.hpp:3008-3011). `bcf_get_format_string` copies the binary GT bytes (src/htslib/vcf.c:3772-3809), and htslib encodes phased ALT as 5 and unphased ALT as 4 (src/htslib/htslib/vcf.h:608-609). As a result, `0/1` is lost entirely and `1/1` counts only on HAP1. Tested: `0/1` gave nothing, `1/1` gave HAP1 only, `0|1` gave HAP2. Keep `-f 1` random phasing (the default, sv_call_main.hpp:133).
4. Two SVs are "near" if any begin/end pair is within `Distance` (analysis.hpp:2845-2853). Every SV with at least one near neighbour is printed as TRUE_POS plus its 500-bp bins (analysis.hpp:2869-2880).
5. Only chromosome indices 0-23 are processed (analysis.hpp:2857, 2974-2975). Index ≥24 is written out of bounds (analysis.hpp:3009).

### 1.2 Contig clustering (`FC_joint_ana`, FC_joint_calling.hpp:2982-3054 → `contig_clustering` 2749-2914)

| Stage | What the code does | Lines |
|-------|--------------------|-------|
| Load | Strict 2-line FASTA (header, sequence). Record 0 must be named `ref*`. The FULL_COVER read number is set to a constant 99. | 794-806, 3033 |
| HOM/HET | Groups contigs by the name prefix before `_chr`. A sample with a single contig is marked HOM and counted as 2 haplotypes. | 786-792, 808-827, 2464 |
| RM parse | Skips 3 header lines; columns 5/6/7/10/11 = query, begin, end, repeat, class/family | 837, 133-147 |
| TRF parse | `@name` blocks; repeats shorter than 30 bp are discarded; ±10 bp edge padding | 853-862, 214-219, 29 |
| Contig→ref | Global ksw2 (`ksw_extd2_sse`, Contig_aligner.hpp:164-177), params `1,4,18,3,65,1` (match, mismatch, open, ext, open2, ext2; Contig_aligner.hpp:91-97), zdrop = max length. Must align end to end. Contigs >30 kb are skipped. | 682-724 |
| VNTR mask per contig | TRF intervals plus RM classes containing `Simple_repeat` or `Satellite` | 624-638, 165-169 |
| Pairwise compare | Greedy: each contig is aligned to the main contig of each existing cluster, in sorted order (cluster 0 = ref). Gaps ≥50 bp inside the ≥50-bp-anchored span become SV_ITEMs. | 2794-2821, 926-938, 1145-1219, 2013 |
| SV typing | 1 = VNTR-ins-in-VNTR, 2 = novel-in-VNTR, 3 = VNTR-ins-in-non-VNTR, 4 = novel-in-non-VNTR. "Novel" means the non-VNTR part is >200 bp **or** >60% (ONT) / >20% (HPRC) of the insert. | 940-944, 1315-1336 |
| MEI exception | In the insert, MEI bases (RM class contains SINE/LINE/Retroposon/DNA/RC/LTR) are counted as non-VNTR unless TRF motifs <6 bp cover ≥80% of them | 172-180, 520-574 |
| Re-typing | Type 3 → 1 if within 30 bp of a ref VNTR with ≥5 shared 10-mers. Type 3 → 4 if covered by a TRF record with period ≥50 and ≤2.5 copies (a DUP). | 1347-1381, 1263-1281 |
| TR vs non-TR | A contig **joins** a cluster only when every SV against its main contig is type 1, i.e. TR length polymorphism is tolerated and any non-TR SV opens a new cluster | 1006-1010, 1384-1385, 2813-2819, 2834 |
| Sub-region pass | Re-clusters within FULL, each ref VNTR region and each ref non-VNTR region. Writes `AT_subregion ... sv_type_count: [t0..t4] ...` lines. | 2291-2318, 2491-2604, 2647-2667 |

### 1.3 Outputs of `FC_joint_ana` (all in CWD)

`CSV.main.log`, `sub_region_detail.log`, `csv_subregion_detail.log`, `cluster.main.log`, `cluster.ANNO.RM.log`, `cluster.ANNO.TRF.log`, `contig_VNTR_AND_ANNO.log`, `var.main_NC.log`, `var.anno_TRF_RM_NC.log`, `contig.main_NC.log`, `var.main.log`, `var.anno_TRF_RM.log`, `contig.main.log`, `nestedStructures.log` (FC_joint_calling.hpp:2862-2913). There is also heavy stdout/stderr logging. Tested: all 14 files appeared.

### 1.4 `contig_clustering_data_preprocessing.py`

- Parses `AT_subregion` lines: ref name, `sv_type_count`, `cluster_size`, `type_VINV_GREEN_sv_base_absolute_BP_in_ref`, `cluster_size_with_ONLY_Insert_VNTR_svs` (lines 87-103). It drops the first AT_subregion line per sample (lines 76-77, 156-157); per FC_joint_calling.hpp:2301-2306 that line is the FULL region.
- Colour codes come from the counts of types 1-4 (lines 6-28, 133-138). In VNTR regions: yellow = type 1 only, green = any type 2. In non-VNTR regions: purple = type 3 only, blue = type 4 only, red = both.
- Corrects cluster_size by merging type-3 breakpoints within 20 bp (`simple_cluster`, lines 31-59, 140-147).
- Extracts the region from the regex `chr[0-9XYM]+_(\d+)_(\d+)_(\d+)+_(\d+)_...`, with a split-on-`_` fallback (lines 110-130). The pattern is human-style.

## 2. MEI-TR array detection (`py/mei_tr_array_single_region_detection.py`)

**Invocation:** `python3 mei_tr_array_single_region_detection.py -f <RM.out> <TRF -ngs file>`. The script requires exactly 3 arguments (lines 400-412). The README form `python3 ... region.fa.out` (README.md:395) only prints usage (tested). `-l` mode is broken: line 391 calls `analyze_mei_periodicity(file)` without `trf_data` (tested: TypeError).

| Item | Value | Lines |
|------|-------|-------|
| Input | RepeatMasker `.out`, first 3 lines skipped; split on whitespace, needs ≥15 columns | 174, 185-190 |
| Fields used | col 5 query (= "sample"), col 6 begin, col 7 end, col 11 class/family. Strand (col 9), repeat name (col 10) and ID (col 15) are not used. | 197-200 |
| MEI classes | substring match on `['SINE', 'LINE', 'Retroposon', 'DNA', 'RC', 'LTR']`; `mei_type` = first keyword that matches | 205-215 |
| Grouping | per query sequence, then per `mei_type`. All SINE families in one contig are pooled. | 218-224, 241-247 |
| Minimum count | ≥3 MEIs per sample (233) and per type (248) | 233, 248 |
| Length | `end - start + 1` | 209 |
| Spacing | `start[i] - start[i-1] - length[i-1]` between consecutive same-type elements; the mean is appended so the lists have equal length | 257-270 |
| DBSCAN | `eps=0.5, min_samples=3` on StandardScaler-standardised (length, spacing) pairs; returns nothing if ≤2 points | 120-143 |
| Cluster pick | largest non-noise cluster | 279-292 |
| Criteria | ≥3 members (295); **reject if** `length_cv > 0.2` **or** (`distance_stdev > 15` **and** `distance_cv > 0.2`) (312). Accepts when length SD ≤20% of mean AND (spacing SD ≤15 bp OR spacing SD ≤20% of mean). Confirmed. | 295, 303-313 |
| Unused | `avg_distance_for_CV = max(avg_distance, 50)` is computed but never used | 306 |
| TRF overlap | intervals from the TRF `-ngs` file that fully contain the array, plus % coverage | 9-51, 53-97, 339-340 |
| Outputs | one stdout line per array (Sample, Main_MEI_type = dominant class/family, cycles, mean length/spacing, CV, SD, TRF intervals, start/end, coverage%) plus a summary (365-378). Writes `mei_data.txt` to CWD (sample, length, class, start, end, spacing) at the very end (381-383). | 363-383 |
| Crash | If an array's sample has no `@` block in the TRF file, `calculate_coverage` returns `[]` (78-80) and `:.2f` raises (376). The run aborts before `mei_data.txt` is written (tested). | 78-80, 376 |

**Changes needed for a custom SINE library or reptile classes**

1. Run RepeatMasker with `-lib your_SINEs.fa` (headers `>name#SINE/xxx`) or a suitable `-species`. Any class string that contains `SINE` is picked up (line 205). RepeatModeler `Unknown` families are ignored.
2. To separate SINE families (e.g. tRNA-derived vs 5S-derived), change the grouping key from `mei_type` (line 243) to `columns[9]` (repeat name) or to the full class/family.
3. Satellite arrays (e.g. CLsat) are **not** MEIs here. Add `'Satellite'` to `mei_keywords` (line 205), but note that `Satellite` also counts as VNTR in the C++ clustering (FC_joint_calling.hpp:165-169).
4. Fragmented RM hits (one SINE split over several `.out` lines) count as separate elements. Merge by column 15 (ID) before running.
5. Make line 340 robust: return 0 instead of `[]` at line 80.
6. Fix the `-l` mode at line 391.

## 3. Contig clustering: who does it

All clustering is in C++: `FC_joint_calling_handler` in src/SVcalling_core/ANALYSIS/FC_joint_calling.hpp, reached through `HitSV tools FC_joint_ana` (analysis.cpp:132-135, 167). It is greedy, reference-seeded and pairwise. ksw2 global alignment is done twice: each contig vs the reference (FC_joint_calling.hpp:684-724, for absolute coordinates and nesting), and each contig vs each cluster's main contig (FC_joint_calling.hpp:926-938, 2809-2812). Clustering only tolerates type 1 (TR-within-TR) differences (see 1.2). The Python preprocessing script only summarises the logs and does no clustering of its own beyond merging breakpoints within 20 bp.

## 4. FC_BED (`-b`) force-contig mode

| Aspect | Finding | Lines |
|--------|---------|-------|
| Option | `-b/--FC_BED`. Help text: "Only support LRS_BAM (-l)", "shorter then 50,000 bp" | sv_call_main.hpp:166-169 |
| Loading | Tab-split, ≥3 fields. **Only contigs with FASTA index 0-23 get coordinates**; other entries stay uninitialised. | SV_core.hpp:862-884 |
| 50-kb limit | **Not checked on the BED.** In PURE_STR, an assembled contig >50 000 bp is skipped (`SKIP: contig length is over long`), but the ASM preset raises the cap to 3 000 000 | SV_core_LRS.cpp:2324-2328, 2367-2371, 2386-2387 |
| LRS/ASM only | FC runs only in the germline LRS branch (`!is_somatic_calling && is_FC_REGION_LRS`), using `LRS_read` signals. The normal pipeline (including SRS) is skipped whenever a BED is loaded. Behaviour with `-n` only and no `-l` is unverified. | SV_core.hpp:2358, 2464-2492 |
| Region placement | A BED interval is used only if it lies **entirely inside** the current 5-Mb reference segment (`SEGMENT_LEN 5000000`). Intervals that straddle a segment boundary are silently ignored. | SV_core.hpp:886-895; RefHandler.hpp:30, 243-252 |
| Signal | Each region becomes one `LRS_SIG_X` signal with support forced to 1000 | SV_core.hpp:889-893 |
| Load margin | 500 bp each side | SV_core.hpp:2230; SV_core_LRS.cpp:2331 |
| PURE_STR output | No VCF header (SV_core.hpp:2314). PASS line: `REF_REGION=chr:st-ed  STAT=PASS  GT:QU:GQ:GI:CI:CN=…  SUP_R:UNSUP_R:FULL_COVER_R:CLIP_AT_LEFT:CLIP_AT_RIGHT:IS_FINAL=…  CONTIG_STR=…`, one line per selected haplotype contig. Failure lines: `STAT=NO_ENOUGH_READ_SUPPORT`, `STAT=ASSEMBLY_FAIL`, `STAT=NO_CONTIG_SELECTED`. | SV_core_LRS.cpp:2338-2407 |
| VCF output | Normal VCF (SV + SNP/INDEL), INFO `SIGNAL_REGION`. `CONTIG_POSITION`, `CIGAR` and `CONTIG` are added only when \|SVLEN\| ≥ 10. Tested on demo. | NovaSVRst.hpp:687-703 |
| `-B` values | Only `PURE_STR` and `VCF` are recognised. Anything else (e.g. the README's `-B BED`, README.md:365) leaves OUT_MODE_UNSET, and the FC branch then writes nothing. | sv_call_main.hpp:237-242; SV_core.hpp:2476-2479 |
| Index build | Ref k-mer index is skipped when not in VCF mode | sv_call_main.hpp:304-305 |

## 5. Dependencies and hardcoded values

| Kind | Item | Where |
|------|------|-------|
| External tools | RepeatMasker (no `-species/-lib`), TRF with `2 7 7 80 10 10 500 -h -d -ngs`, bedtools, SLURM `sbatch` | cohort_csv_analysis.md:144-145, 228-239, 86 |
| Python | numpy, scikit-learn (DBSCAN, StandardScaler) | mei_tr_array_single_region_detection.py:5-7 |
| Python | pandas, matplotlib, gzip | mei_tr_array_gene_regulatory_annotation.py:2-5 |
| Python | stdlib only (re, os, argparse) | contig_clustering_data_preprocessing.py:2-4 |
| C++ | `<filesystem>` (C++17) | FC_joint_calling.hpp:13, 22 |
| Reference | `GRCh38.fa` in every doc script | cohort_csv_analysis.md:71, 104, 162, 227 |
| Annotation | `Homo_sapiens.GRCh38.regulatory_features.v115.gff3.gz`, `gencode.v44.basic.annotation.sorted.gtf.gz` (CWD) | mei_tr_array_gene_regulatory_annotation.py:238-239 |
| Fixed output names | `mei_data.txt`; 14 FC_joint_ana logs in `./` | mei_tr_array_single_region_detection.py:381; FC_joint_calling.hpp:2997 |
| Developer paths | `/media/fenghe/...` in unregistered `vcf_add_alt_string` | analysis.hpp:2241-2244, 2320 |
| Developer path | `/home/fenghe/E` in `forceCalling_run`, reachable via `HitSV ngs_fc` | forceCalling.hpp:860, 873; main.cpp:24-26, 113 |
| 24-chromosome caps | csv_region_500w, FC_BED loader, entropy/small_var_den loaders, genomic_contest_* | analysis.hpp:2857, 2974; SV_core.hpp:878; analysis.hpp:65, 538; Genome_context_analysis.hpp:220, 251, 268, 986 |
| Human TE names | None in the population/analysis code. MEI is RepeatMasker **class** substrings only (above). The only Alu-specific path is SRS `force_calling_ALU` (on by default) | sv_call_main.hpp:175; SV_core_SRS_ins_del.cpp:2118-2120 |
{% endraw %}
