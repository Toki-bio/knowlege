---
title: "hitsv: Control flow: per-block pipeline and LRS mode"
type: reference
generated: true
---

> Part of the [hitsv](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/hitsv/references/control-flow-lrs.md)

{% raw %}
# Control flow: per-block pipeline and LRS mode

HitSV v2.0.4, commit 93880367ef76f6440c918f7dd9744dda2659bcdb. All paths relative to `src/`. Every line number below was checked by reading the source.

## 0. Block loop (context)

| Step | Location | What happens |
|---|---|---|
| Block size | `SVcalling_core/RefHandler.hpp:30` | `SEGMENT_LEN 5000000`: blocks are 5 Mb, not 2 Mb (the comment at `sv_call_main.hpp:298` says "2M"). |
| Next block | `RefHandler.hpp:243-267` | Blocks do not overlap (`st_pos + SEGMENT_LEN`, :248). The whole chromosome is loaded into memory once per chromosome (`load_reference`, :306-327). |
| Ref k-mer index | `sv_call_main.hpp:304-305` -> `RefHandler.hpp:273-291` | Built only in VCF output mode. k=16 (`RefHandler.hpp:209`), sampling every 3rd position (:283). Lookups return a hit only for k-mers seen exactly once (:297). |
| Pipeline call | `sv_call_main.hpp:308` | `SV_calling_in_region_main_pipeline()` |

## 1. `SV_calling_in_region_main_pipeline` (`SVcalling_core/SV_core.hpp:2355-2637`)

### Decision points
1. `is_somatic_calling = (LRS_Tumor_BAM_F != NULL)` (`SV_core.hpp:2356`). The tumor BAM comes from `-X` (`sv_call_main.hpp:112`).
2. `is_FC_REGION_LRS` is set when a `-b` BED was loaded (`SV_core.hpp:2358`, `sv_call_main.hpp:169,287-288`).
3. Otherwise the basic path runs: LRS if `-l` is given (`SV_core.hpp:2495`), then SRS if `-n` is given (`SV_core.hpp:2567`). Hybrid means both blocks run.

### 1a. Somatic branch (`SV_core.hpp:2360-2460`)
1. Collect tumor signals with `is_LRS_somatic=true`, then combine them (`:2369-2370`).
2. Minimum tumor support is `depth*0.1`. There is no `MAX(...,1)` floor here (`:2371`).
3. Collect and combine normal signals (`:2374-2375`). The normal BAM is optional (`-x`, `sv_call_main.hpp:113`). If no normal BAM is given, normal signal collection still runs on the uninitialised normal handler; this path is unverified.
4. Maximum reads for assembly: tumor `MAX(10*depth,100)`, normal `MAX(8*depth,100)` (`:2378-2379`).
5. For each tumor signal that overlaps the block and has enough support (`:2387-2394`):
   - assemble the tumor reads (`:2399`);
   - assemble the normal reads (`:2410`). **If normal assembly fails, the signal is skipped** (`:2417`);
   - call `LRS_SV_generating_Somatic` (`:2426`).
6. Then: `LRS_remove_duplications`, `LRS_combine_Homozygous_SVs`, `calGT_LRS`, optional random phasing, and VCF output (`:2435-2459`). SRS re-genotyping is not done (`:2449`).

### 1b. Force-calling (BED) branch (`SV_core.hpp:2464-2490`)
1. Collect signals (`:2471`).
2. **Signal combining is skipped.** The BED regions are injected directly as signals (`:2474`, `SV_core.hpp:886`).
3. Run `LRS_SV_CALLING_germline` (`:2477`). With `-B PURE_STR`, contigs are printed instead of VCF records (`SV_core_LRS.cpp:2354-2407`).

### 1c. Basic germline branch (`SV_core.hpp:2492-2636`)

**LRS-only or the hybrid LRS part (`:2495-2565`):**
1. `LRS_SIG_COLLECTION` (`:2502`) -> `LRS_SV_SIG_COMBINE` (`:2505`) -> `LRS_SV_CALLING_germline` (`:2506`).
2. `LRS_remove_duplications` (`:2511`) and `LRS_combine_Homozygous_SVs` (`:2515`). Two 0/1 calls of the same variant on different haplotypes of one region become 1/1.
3. `calGT_LRS` (`:2519-2520`).
4. **Hybrid only:** re-genotype each LRS SV with SRS reads (`:2523-2549`), subject to these rules:
   - skip if `depth<10`, preset is ERR_PRONE, and the site is a VNTR according to the `-I` local-repeat index (`:2529-2536`);
   - skip if `|SVLEN| < MIN_sv_len` (`:2537`);
   - force 0/0 if `supp==1 && not_supp>=1` and SRS REF reads >= 20x ALT (`:2545-2546`; `NovaSVRst.hpp:604-613`);
   - replace GT with the SRS GT if `not_supp==0 && supp<4` (`:2547-2548`; `NovaSVRst.hpp:599-601`).
5. Random phasing (on by default, `-f 1`, `sv_call_main.hpp:133`) and VCF write (`:2554-2564`). Every VCF line is also echoed to stderr (`:2562`).

**SRS-only or the hybrid SRS part (`:2567-2635`):**
1. With `-Y` (INNER_TEST) and LRS present, SRS calls are loaded from a VCF in debug mode (`:2579-2583`).
2. Otherwise the full SRS pipeline runs (`:2588-2608`): signals, clustering, small, VNTR, large-DEL, dedup, GT, INV, BND.
3. **Hybrid only:** `Hybrid_regenotype_NGS_USING_TGS_data` (`:2612-2615`, `SV_core_LRS.cpp:2636-2687`):
   - an SRS INDEL within 300 bp of any LRS SV is dropped (`SV_core_LRS.cpp:2616`);
   - the remaining SRS INDELs are re-genotyped by counting LRS reads that fully span ±500 bp. More than 3 such reads gives 0/0; 1 to 3 gives 0/1 (`:2659-2677`). This count is coverage, not support for the variant.
4. Write the SRS INDEL, INV and BND results (`:2626-2633`).

| Mode | Steps executed |
|---|---|
| LRS-only (`-l`) | 1c LRS part; SRS part skipped |
| SRS-only (`-n`, needs `-I`, `sv_call_main.hpp:197-200`) | 1c SRS part only |
| Hybrid (`-l -n`) | LRS part + SRS re-GT of LRS calls, then SRS part + SRS-vs-LRS dedup/re-GT. A warning is printed if `-k` is missing (`sv_call_main.hpp:259-261`) |
| Somatic (`-X [-x]`) | 1a only |
| BED force (`-b`) | 1b only |

## 2. LRS signal collection: `LRS_SIG_COLLECTION` (`SVcalling_core/SV_core_LRS.cpp:1945-2178`)

### Preset switches (`:1951-1962`)
| Flag | ONT_Q20 | HIFI | ERR_PRONE | ASM |
|---|---|---|---|---|
| dense small-variant signal (`using_small_signal`, :1954-1955) | on | on | off | off |
| clip signal (`using_clip_signal`, :1958-1960) | on | on | off | on |
| SA parsing (`using_SA_signal`, :1961-1962) | on | on | off | on |
| SNP signal | always off (`:1956-1957`, "not used now") | | | |

### Read filters
- Secondary alignments are skipped (`:1994-1995`); MAPQ < 1 is skipped (`:1996`); supplementary alignments with `l_qseq < 1200` are skipped (`:1999`).
- Read QC applies to every preset except ASM (`:2004`; `LRS_read_quality_check` at `:2468-2528`):
  - a read fails if any 1000-match window has error > 0.20 (`:2498-2505`);
  - a read fails if its overall error > 0.04, or > 0.16 under ERR_PRONE (`:2519-2522`).
  - Error here counts mismatches plus one per I/D op. The `-e` option is parsed (`sv_call_main.hpp:146`) but never used (grep finds only the declaration and the parser).

### Signatures and thresholds
| Signature | Rule | Line(s) |
|---|---|---|
| CIGAR INS | `len >= MIN_sv_len` -> `LRS_SIG_INS` at the ref position | `:2066-2076` |
| CIGAR DEL | `len >= MIN_sv_len` -> `LRS_SIG_DEL` over [pos, pos+len] | `:2082-2093` |
| `MIN_sv_len` | `-m`, **default 25** | `sv_call_main.hpp:134`; `SV_core.hpp:2300` |
| Soft clip | `length > 1000` -> CLIP_LEFT/RIGHT signal | `:2038-2064` |
| Hard clip | `> 1000` sets `with_long_CLIP`; it is stored as a signal **only in somatic mode** | `:2114-2132` |
| Dense small variants | per read, sliding 200 bp window over sub-threshold I/D (weight = length) and mismatches (weight 1); emit `LRS_SIG_X` when the queue has **>10 events or >35 weighted bp** | `:2077-2080, 2094-2097, 2103-2105`; `:2535-2555` (window :2539, trigger :2547) |
| SA / split alignment | parsed only for reads with a >1000 bp clip (`:2144-2160`). The `sa_handler` is read only in somatic BND generation (`SV_core_SRS_ins_del.cpp:3758`). **In germline mode SA tags produce no signal.** | |

### Depth estimate (`:2166-2177`)
`average_read_depth = total_primary_bases / (ave_read_len + block_len)`. It is computed over the 5 Mb block, so it is a block-level mean, not a local depth.

### Local repetitiveness R
In LRS-only runs, `rlr` stays NULL because it is allocated only when `-n` is given (`sv_call_main.hpp:256,263-271`). In the LRS path its only use is the hybrid ERR_PRONE re-GT skip (`SV_core.hpp:2529-2535`). **R plays no role in LRS signal detection, assembly or calling;** it is an SRS and hybrid feature.

## 3. Candidate region formation: `LRS_SV_SIG_COMBINE` (`SV_core_LRS.cpp:2180-2277`)

| Step | Rule | Line(s) |
|---|---|---|
| Support floor | `min_support = MAX(int(depth*0.1), 1)` | `:2189-2190` |
| P0: small-signal cluster | sorted; consecutive signals merge if both POS and END are within 20 bp (`similar_SIG_SMALL`); the chain breaks at the first non-match; keep if `support >= min_support && >= 2` | `:2194, 2197-2218`; `SV_core.hpp:447-455` |
| P1: CIGAR+CLIP cluster | CIGAR signals appended to CLIP; consecutive signals merge if any endpoint pair is <300 bp apart, or one contains the other (`similar_SIG_Normal`). **Signal type is not compared**, so INS, DEL and CLIP merge together. Keep if `support >= min_support` | `:2192, 2221-2246`; `SV_core.hpp:432-445` |
| P2: region merge | P1 output plus P0 output; merge if the regions overlap after ±600 bp padding (`region_overlap`); support counts are **summed**; inner loop has no break (O(n²)) | `:2193, 2248-2269`; `SV_core.hpp:457-467` |
| Germline gate | `support >= int(depth*0.1)` (no floor of 1 here) | `:2316, 2338` |
| Out-of-block | skip if the signal does not overlap the current block | `:2334-2336` |

**Maximum region length:** the code has no explicit cap on the merged signal region. Merging is transitive (each `l[i]` keeps growing), so dense repeat-rich loci can produce very long regions. The practical limits are downstream:
- read-count cap `MIN(MAX(5*depth,100),1000)` (`:2320-2321`), enforced at `SV_core_SRS_ins_del.cpp:2969-2973`;
- 50 kb contig cap in PURE_STR output only (`:2324, 2367`); see pitfalls.

**Extension around the signal** (`LRS_assembly_variations`, `SV_core_SRS_ins_del.cpp:2854-2985`):
- `edge_len = sig_para.MaxReadLen` (`:2863`) is subtracted **from the left only** (`:2872-2879`). `sig_para` is zeroed (`SV_core.hpp:1929`) and filled only when an SRS BAM is present (`SV_core.hpp:2262-2264`). So edge_len is 0 in LRS-only runs and equals the SRS read length in hybrid runs.
- Reference window = [POS - edge_len - **500**, END + **500**] (`refAdditionLoad=500`, `:2864, 2903-2906`).
- Read window = same ±**500** (`readAdditionLoad=500`, `:2865, 2952-2958`).
- Minimum read slice length `readAdditionLoad*2-200 = 800` bp. The comment at `:2961` says 1800, which is wrong.
- Reads clipped >200 bp at both ends are excluded (`SV_core_SRS_ins_del.cpp:1013-1014`); supplementary alignments with MAPQ < 20 are excluded (`:1005-1006`).
{% endraw %}
