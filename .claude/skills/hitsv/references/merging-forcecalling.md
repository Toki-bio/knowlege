# ngs_merging and ngs_fc

HitSV v2.0.4 (commit 93880367). Both are registered in `src/main.cpp` with the same help text: "Force calling tools for NGS dataset" (main.cpp:112-113). Dispatch:
- `ngs_merging` -> `NGS_SV_MERGING_HANDLER::SV_merging` (main.cpp:30-34)
- `ngs_fc` -> `FC_handler::forceCalling_run` (main.cpp:24-28)

**Bottom line: both are author-internal 1KGP research scripts with hardcoded absolute paths. They do not work as general tools.**

## ngs_merging (`src/SVcalling_core/NGS_merging.hpp`)

| Aspect | What the code does | Lines |
|---|---|---|
| CLI arguments | `argc/argv` are accepted but **never read** | NGS_merging.hpp:182-185 |
| Reference | Hardcoded `/media/fenghe/Data/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa` | NGS_merging.hpp:184 |
| Sample list | Read from the `#CHROM` header columns 9+ of hardcoded `/media/fenghe/MyPassport_4T/bamdata/1KGP_HitSV/filltag_S5.vcf` (only the first 1000 lines are loaded) | NGS_merging.hpp:185, 191-210 |
| Per-sample input | Hardcoded `/media/fenghe/MyPassport_4T/bamdata/1KGP_HitSV/sorted_vcf/<sample>_sv_spa.sort.vcf.gz` | NGS_merging.hpp:220-227 |
| Contig limit | Index built for at most 24 contigs; records with chrID>=24 are skipped (human chr1-22,X,Y) | NGS_merging.hpp:171, 239-240 |
| Record filter | Skip `BND`; keep FILTER id 0 or 13 only | NGS_merging.hpp:241-245 |
| Bucketing | 30-bp position blocks (`SV_MERGE_IDX_B_SIZE = 30`). **`init()` is never called by `SV_merging`**, so the value is uninitialised when used (no call to `init()` found in the file) | NGS_merging.hpp:163-166, 178, 258 |
| Match rule (`sameWith`) | Requires an **exact same POS** in the same 30-bp bucket. Then: ALT and alt <400 bp needs an exact SVLEN match; otherwise SVLEN ratio in (0.9,1.1). If the bucket holds >100 SVs ("high complex"), the rule switches to ALT-length ratio in (0.9,1.1). | NGS_merging.hpp:98-109, 258-263 |
| Output (stdout) | Text, not VCF: `chrID:pos  SVLEN  REF  ALT  SA:` followed by tab-separated `sampleID:GT1:GT2:QUAL` (raw htslib-encoded GT values) | NGS_merging.hpp:284-300 |
| Purpose | Population-level union/merge of per-sample HitSV SRS VCFs into a site list with carrier samples | inferred from NGS_merging.hpp:212-303 |

Unused code: `sameWithM2` (NGS_merging.hpp:111-124) and `count_sample_region` (79-95). Neither is called within the file.

## ngs_fc (`src/SVcalling_core/forceCalling.hpp`, `forceCalling.cpp`)

| Aspect | What the code does | Lines |
|---|---|---|
| Positional args | `argv[1]` = ref.fa, `argv[2]` = VCF list, `argv[3]` = BAM, `argv[4]` = local-repeat file, `argv[5]` = sampleID | forceCalling.hpp:863-867 |
| Arguments actually used | Only ref (888), BAM (883) and sampleID (925). `vcf_fn_in` and `local_repeat_fn` are read but never used. | forceCalling.hpp:864, 866 |
| SV site input | Hardcoded path `/home/fenghe/E`, up to 50000 lines, in the **`ngs_merging` output format** (`chrID:pos SVLEN REF ALT SA: id:GT1:GT2:QUAL ...`) | forceCalling.hpp:873-880, 819-834 |
| Hard position skip | `if(SV_info[SV_ID].pos < 23504400) continue;` is a leftover debug coordinate | forceCalling.hpp:935 |
| Clustering | Sites within 30 bp form one cluster. Classes: SVLEN>800 long INS; SVLEN<-700 long DEL; REF and ALT both >30 complex; otherwise small | forceCalling.hpp:903-922 |
| Force-call condition | The target sample is **not** already a carrier, and the cluster has more than one carrier entry | forceCalling.hpp:924-931, 939 |
| Genotyping | Small clusters only (1-19 haplotypes): REF+ALT haplotypes are built over +/-300 bp of reference (892, 942-966). Reads are realigned to the haplotypes with `GT_LOCAL_handler` (979; 308-317). | forceCalling.hpp:892, 942-979 |
| GT rule | 1/1 if sup > 4*ref; 0/0 if sup*3 < ref; else 0/1. `./.` if unknown > 3*(sup+ref), or if both sup<4 and ref<4. | forceCalling.hpp:371-380 |
| Long INS/DEL | Long INS/DEL and complex clusters are collected but never genotyped. The DRP genotyping block is commented out. | forceCalling.hpp:983-997 |
| DRP helpers | Hardcode `max_isize = 1000` and `normal_read_length = 150` | forceCalling.hpp:709-710 |
| Read length | `gth.init(&c_b, 150, ...)` uses a hardcoded 150 | forceCalling.hpp:889 |
| Output | **stderr log only**. `GT_main` returns a GT that `forceCalling_run` discards (979), and nothing is written to stdout or to a VCF. `showGT` prints to stderr. | forceCalling.hpp:308-316, 979; forceCalling.cpp:13-22 |
| Purpose | Population force-calling (re-genotyping other samples' SV alleles in one sample's SRS BAM), the 1KGP second pass | inferred from forceCalling.hpp:860-1005 |

## Practical use for non-human projects

Neither subcommand is usable without editing the source:
- **Hardcoded paths:** NGS_merging.hpp:184-185, 223; forceCalling.hpp:873.
- **Human 24-contig cap:** NGS_merging.hpp:171, 239.
- **Leftover debug coordinate:** forceCalling.hpp:935.
- **No output:** `ngs_fc` produces no results file (forceCalling.hpp:979).

For multi-sample merging, use an external merger (e.g. Jasmine, SURVIVOR, truvari collapse) on per-sample `call` VCFs. This is a recommendation, not something the code provides.
