# `HitSV tools` sub-tools

The `tools` command is registered at src/main.cpp:115 (`ch.add_function("tools", "Small tools", analysis_main)`). `analysis_main` builds its own dispatcher at src/SVcalling_core/ANALYSIS/analysis.cpp:156-196. That dispatcher calls `set_main_command("analysis")` (analysis.cpp:159), which is why `HitSV tools` prints "Usage: HitSV analysis <command>". Commands match exactly and are case-sensitive (src/cpp_lib/get_option_cpp.hpp:367). The built binary also prints `Version: 2.0.0` (src/clib/desc.h:2) even though the release is v2.0.4.

All file paths below are relative to `src/SVcalling_core/ANALYSIS/` unless shown otherwise.

## Registered (active) sub-tools: 13

| Sub-tool | Registered | Implementation | Arguments (real) | What it does / output |
|----------|-----------|----------------|------------------|-----------------------|
| `csv_region_500w` | analysis.cpp:165 | analysis.hpp:2947-3018 | `VCF Distance min_SV_len` | Marks PASS SVs on HAP1/HAP2 that lie within `Distance` of another SV, then prints `chr st ed HAPn TRUE_POS` plus 500-bp bins `… 500W` to stdout (analysis.hpp:2873-2879). Drops DELs because SVLEN<0 (analysis.hpp:2999), uses phased GT only (3008-3011), and handles chr index ≤23 only (2857). |
| `contig_file_split` | analysis.cpp:166 | analysis.hpp:3050-3137 | `list_of_PURE_STR_files work_dir ref.fa` | Writes one FASTA per region (`<chr>_<st>_<ed>`) into work_dir: the ref sequence first, then `<sample>_<region>_<CI>` contigs. The sample name is taken from the parent directory (3077-3078). |
| `FC_joint_ana` | analysis.cpp:167 | FC_joint_calling.hpp:2982-3054 | `region.fa RM.out TRF.ngs.out HPRC\|ONT` | Population haplotype clustering with RepeatMasker/TRF-aware SV typing. Writes 14 logs to CWD, including `CSV.main.log` and `nestedStructures.log` (2862-2913). The help string "[TRUR/FALSE]" is wrong: the code requires HPRC/ONT (3003-3014). |
| `dump_ref_by_region` | analysis.cpp:169 | analysis.hpp:1545-1585 | `ref.fa region_file` (tab chr st ed; analysis.hpp:1494-1527) | Prints FASTA `>chr_st_ed len_…` for each region, 70 columns per line (1572-1579). |
| `get_ref_region` | analysis.cpp:170 | analysis.hpp:1408-1421 | `ref.fa chr:st-ed` | Prints the sequence via faidx (tested). |
| `ref_split` | analysis.cpp:171 | analysis.hpp:1435-1465 | `input.fa` | Writes `<seqname>.fa` per sequence and prints the headers to stdout (1449-1454). |
| `getReverseStr` | analysis.cpp:172 | analysis.hpp:1599-1609 | `ACGT…` | Reverse complement, printed to **stderr** (1607; tested). |
| `entropy_analysis` | analysis.cpp:173 | analysis.hpp:469-515 | `region_bed ref.fa vcf bam min_sv_len` (help 455). The registered help "[bam_file] [ref.fa]" is wrong. | Per BED region, uses LRS breakpoint signals in 10-bp windows plus VCF SVs and prints `region total_KLD n_windows` (391-451). The BED loader skips line 1 (61). |
| `small_var_den` | analysis.cpp:174 | analysis.hpp:694-754 | `region_bed ref.fa bam` (help 682) | Per region, prints identity and mismatch/indel density lines `R_FULL` and `R_PART` (±500 bp trimmed) (724-747). The BED loader skips line 1 (534). |
| `align_seq` | analysis.cpp:176 | FC_joint_calling.hpp:2925-2979 | `target query match mismatch go ge go2 ge2 FULL\|ALN\|CIGAR` | ksw2 global alignment of two strings using the same aligner as FC_joint_ana. The argc check is `<8` but the code reads argv[9] (2927, 2942), so it segfaults when the mode argument is missing (tested). |
| `genomic_contest_Analysis_SRS` | analysis.cpp:177 | Genome_context_analysis.hpp:559-~790 | `ref.fa bam block_len step` (exactly 5 argv; 561) | Per-window SRS signal counts (CIGAR/CLIP/DRP/small-var). Output table header at Genome_context_analysis.hpp:746. |
| `genomic_contest_Analysis_LRS` | analysis.cpp:178 | Genome_context_analysis.hpp:793-948 | same | Same table for LRS (943 → 739-760). The window end is hardcoded as `pos+600` whatever the block length (755). |
| `genomic_contest_Analysis_ANNO` | analysis.cpp:179 | Genome_context_analysis.hpp:950-1007 | `ref.fa vcf(s) VNTR_file HIGH_CONF.bed` | Annotates windows with SV type, max VNTR value and high-confidence flag. Header at 981; VNTR file format at 236-272 (first line `block_len step`, then `chrID pos maxVNTR`). The help example has an extra `input.bam` (410). |

All `genomic_contest_*` tools register with an empty help string (analysis.cpp:177-179) and cap to 24 sequences (Genome_context_analysis.hpp:220, 251, 986).

## Present in source but not registered (commented out)

| Function | Commented at | Implementation |
|----------|--------------|----------------|
| `vntr_analysis` | analysis.cpp:161 | VNTR_analysis.hpp:398 |
| `bed_in_repeat_masker_region` | analysis.cpp:162 | VNTR_analysis.hpp:154 |
| `csv_region` (non-500W twin) | analysis.cpp:164 | analysis.hpp:2886-2945 |
| `genomic_contest_heatMapData`, `…Sum` | analysis.cpp:181-182 | Genome_context_analysis.hpp:1552, 1774 |
| `fa_stat_full`, `string_local_repeat` | analysis.cpp:184-185 | analysis.hpp:1279, 1356 |
| `vcf_2_contig` | analysis.cpp:187 | analysis.hpp:3139-3224 |
| `bed_merging` | analysis.cpp:190 | analysis.hpp:1176 |
| `T2T_calling_ana`, `t2t_split_alignment_analysis` | analysis.cpp:193-194 | ASM_analysis_handler.hpp:268, 638 |

`SV_in_simple_repeat` (VNTR_analysis.hpp:213), `combine_sort_vcf` (analysis.hpp:774) and `vcf_dump`/`vcf_add_alt_string` (analysis.hpp:2205, 2241) have no wrapper in analysis.cpp at all.

Tested: running `HitSV tools` with no arguments lists exactly the 13 active tools above.
