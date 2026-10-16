# Local assembly, contig realignment, polishing, genotyping

HitSV v2.0.4, commit 93880367ef76f6440c918f7dd9744dda2659bcdb. Paths are relative to `src/`. Every line number below was checked by reading the source.

## 4. Local assembly (LRS): `LRS_assembly_variations` (`SVcalling_core/SV_core_SRS_ins_del.cpp:2854-3218`)

### Which assembler is used
| Assembler | Used by | Evidence |
|---|---|---|
| `MainDBGHandler1` (DBGAssembler.hpp/.cpp) plus read-threading SDP | **every LRS path** (germline, somatic, BED force) | `SV_core_SRS_ins_del.cpp:3010-3013` |
| `MainDBGHandler1` with k=60..120 step 20 (and 65 for error correction) | SRS INS/DEL/VNTR | `SV_core_SRS_ins_del.cpp:2258-2261, 2270-2283` |
| `MainAssemblyHandler` (GreedyAssembler) | **SRS only**: MEI/long SV (`:2025`), INV/TRA (`SV_core_SRS_inv_tra.cpp:416, 674`) | grep finds only these `run_assembly(am)` callers |

### Trace
1. **Load reference and reads.** Windows are ±500 bp (see control-flow §3). Reads are capped at `max_read_depth_assembly`; the region is skipped if more reads load (`:2969-2973`).
   - ASM preset: per read name, the longest same-strand span over primary and supplementary pieces is stitched from the primary sequence (≤3 Mb, `:1029, 1074, 1108`).
   - ASM preset: **only one** sequence is kept, the one farthest from the region midpoint (`:1253-1274`). ASM mode is therefore effectively haploid.
2. **k-mer size** (`:3001-3003`): **k=20**; **k=10 for ERR_PRONE**.
3. **k-mer counting** (`:3010`): only reads are counted, because `using_reference=false` skips mode 0 (`DBGAssembler.hpp:3178-3179`). The reference is not part of the graph.
   - Read-begin and read-end k-mers are recorded only for reads longer than 5000 bp (`DBGAssembler.hpp:3212-3215`).
4. **Low-depth threshold** `minDepth_LRS` (`:2990-3000`). N is the number of loaded reads.
   | Preset | minDepth_LRS |
   |---|---|
   | ONT_Q20 | N/6; if N<10 then (N+2)/3 |
   | ERR_PRONE | N/4 |
   | HIFI, ASM | N/6 |

   A k-mer is dropped when `read_Lrs_n < minDepth_LRS` (and it has no SRS or ref count) (`DBGAssembler.hpp:56-58`; `DBGAssembler.cpp:635-639`). The code has no explicit "high-depth node" threshold.
5. **Unitig dBG** (`DBGAssembler.cpp:623-841`):
   - edges come from (k-1)-overlaps (`:650-681`);
   - unitig breaks at in>1 or out>1 and at read-begin/end k-mers (`:668-708`);
   - unitigs are walked (`:726-744`).
6. **Cycle handling:** a Tarjan-style SCC search marks unitigs in cycles as `isRepeat` (`DBGAssembler.cpp:561-621, 750-775`). **The LRS path never reads `isRepeat`** (no grep hit in `SV_core_LRS.cpp` or `SV_core_SRS_ins_del.cpp`). Repeats are handled by threading reads instead:
   - each read becomes a unitig path (`DBGAssembler.hpp:2710-2800`);
   - same-unitig pieces merge only if the gap offset is ≤50 (`:2440-2442`);
   - in SDP, unitigs occurring >30 times in a target path are ignored as anchors (`SV_core_LRS.cpp:630`).
7. **Read clustering by SDP** (`LRS_contig_candidate_generate`, `SV_core_SRS_ins_del.cpp:2336-2454`):
   - reads are processed full-cover first, then by matched k-mer count (`DBGAssembler.hpp:2998-3003`);
   - reads clipped at both ends are skipped (`:2366-2367`);
   - each read is SDP-aligned to every existing target (`:2399`). Score = matched k-mers − (diagonal shift)², with quadratic begin/end penalties for full-cover reads, and a look-back of 200 matches (`SV_core_LRS.cpp:645-744`, :671, :701);
   - the read is absorbed if `max_score >= MIN_score_final` (**100**, or **−300** for ERR_PRONE) (`:2341-2344`; `SV_core_LRS.cpp:1033`) and no true branch is found;
   - a true branch is: length difference ×5 ≥ max_len and >20 k-mers, or differing unitigs of ≥3 k-mers that are not homopolymer run-length-equal, with total nodes ≤7 and no short (<5 k-mer) nodes (`SV_core_LRS.cpp:802-900`);
   - reads with no target become a new target (`:2411-2417`). Clip mode must match: right clip, left clip, or none (`:2393`).
   - A second pass re-runs this when ≥ `max(minDepth_LRS,2)` reads hit the same "long branch" (>7 nodes) (`SV_core_SRS_ins_del.cpp:2347-2355, 3057`; `SV_core_LRS.cpp:12-40`).
8. **Maximum targets:** germline 1,000,000. **Somatic 15**; exceeding it clears all targets, so the region yields nothing (`:3059-3062, 2418-2420`).
9. **Contig string:** built from the unitig path (`SV_core_LRS.cpp:42-185`). For ERR_PRONE, bases from depth-1 unitigs or read-only gaps get quality '0' (`:134-140, 151-155`). Those bases are N-masked before realignment (`SV_core_SRS_ins_del.cpp:3368-3375`).
10. **Clip contig combination** (germline only): runs when ≥1 left-clip and ≥1 right-clip contig exist, or ≥2 contigs with ≥3 reads including a clip contig (`:3089-3144`; `SV_core.hpp:1313-1317`).
11. **Contig selection** (`CONTIG_SELECTER::select_contig`, `SV_core_LRS.cpp:1208-1300`):
    - support = best-target reads that are full-cover or clip (`:1176-1195`);
    - support is zeroed when average MAPQ < 6 or there is no full-cover read (`:1215-1223`). **This also zeroes pure clip contigs.**
    - contigs whose run-length-compressed lengths are within 5 bp (HIFI) or 10 bp (others) collapse to the best-supported one (`:1230-1254`);
    - a contig starts as final when support ≥3 (`SV_core.hpp:1065`);
    - **if more than 2 final paths remain, the top 2 by support are kept** (`SV_core_LRS.cpp:1265-1269`). This is a diploid cap.
12. **Homopolymer correction** (ONT_Q20 and HIFI only), with k=40 despite the "15" comment (`SV_core_SRS_ins_del.cpp:3170-3188`):
    - unique contig 40-mers bound "bubbles";
    - if ≥2 reads disagree with the contig distance, the median read's sequence replaces the bubble, but only when it differs from the reference distance by ≤10 bp (`:2581-2661`);
    - **the replacement is applied only when exactly one bubble exists** (`:2807-2815`).

### Preset summary (assembly)
| Item | ONT_Q20 | HIFI (`PRESET_CCS`) | ERR_PRONE | ASM |
|---|---|---|---|---|
| k | 20 | 20 | 10 | 20 |
| minDepth_LRS | N/6 (N<10: (N+2)/3) | N/6 | N/4 | N/6 |
| SDP min score | 100 | 100 | −300 | 100 |
| read QC | ≤4% | ≤4% | ≤16% | none |
| clip and small signals | yes, yes | yes, yes | no, no | yes, no |
| contig merge tolerance | 10 | 5 | 10 | 10 |
| homopolymer fix | yes | yes | no | no |
| depth-1 masking | no | no | yes | no |
| small-var output | per `-v` | per `-v` | forced off (`sv_call_main.hpp:219-220`) | per `-v` |
| read loading | slices | slices | slices | 1 stitched seq |

## 5. Contig-to-reference realignment (`SV_core_SRS_ins_del.cpp:3872-4246`; `Contig_aligner.hpp`)

### ksw2 presets (`Contig_aligner.hpp:799-818`)
`set_option(match, mismatch, gapO1, gapE1, gapO2, gapE2, zdrop, band, flag)`

| `-P` | A | B | O1,O2 | E1,E2 | Line |
|---|---|---|---|---|---|
| asm5 (default, `sv_call_main.hpp:153`) | 1 | 19 | 39,81 | 3,1 | `:801` |
| asm10 | 1 | 9 | 16,41 | 2,1 | `:804` |
| asm20 | 1 | 4 | 6,26 | 2,1 | `:807` |
| m2 (hidden; also used for BND aligners `SV_core.hpp:2322, 2325`) | 1 | 15 | 18,65 | 3,1 | `:810` |
| somatic T-vs-N (internal) | 2 | 4 | 4,24 | 2,1 | `:813` |

- Flag is `KSW_EZ_EQX`. The initial z-drop is `gap_open2_D + 600`, but it reads the member before it is assigned (`:801`). Band is 600.
- The aligner is `ksw_extd2_sse` with dual affine gaps (`Contig_aligner.hpp:176-177`).
- Overrides at runtime:
  - `set_region_addition_load_short()` sets `ca.setZdrop(1200,1200)` per signal (`SV_core.hpp:2230`; called at `SV_core_LRS.cpp:2331`).
  - `LRS_assembly_variations` then sets **z-drop 5000, band 5000** (`SV_core_SRS_ins_del.cpp:2921`). `setZdrop`'s second argument is the bandwidth (`Contig_aligner.hpp:565-568`).
- **Long mode** (`Contig_aligner.hpp:181-217`):
  - triggers when |ref−contig| > 1.5×min, or both are >10 kb (`SV_core_SRS_ins_del.cpp:3890-3896`), or when normal mode leaves >500 bp unaligned (`:3916-3923`);
  - pads both sequences with ≤10 kb of the same random sequence on each side (addition_load ≤ `ALN_ADDITION_LOAD_MAX`=20000);
  - band and z-drop are set to **30000** (`:205`), which persists for later calls.
- Long-range contigs (clip-dominated): a 16-mer reference index finds clusters (>60 hits, diagonal within 100) (`SV_core_LRS.cpp:1809-1830`).
  - An SV of 1 kb to 50 kb between clusters (`SV_core.hpp:1606-1608`; `SV_core_LRS.cpp:1671`) defines a two-part reference with 800 bp margins (`SV_core_LRS.cpp:1724-1725`).
  - The middle deletion is then widened by the gap between parts (`SV_core_SRS_ins_del.cpp:3940-3982`).
  - The region is rejected if the middle deletion is >100 kb and either flank has <1500 matched bases (`:3971-3972`).

### Variant extraction (`Contig_aligner.hpp:495-555`, then `SV_core_SRS_ins_del.cpp:3990-4075`)
- `adjustCIGAR` trims small edge matches, ≤15 bp tail (`Contig_aligner.hpp:270-275`).
- Each I or D op with `len >= MIN_sv_len` counts as a true SV. With `-v` (default 1, `sv_call_main.hpp:135`), every I/D and every mismatch is also emitted (`Contig_aligner.hpp:515-523, 530-546`).
- Edge exclusions: first and last CIGAR ops are excluded, and insertions within 20 bp of the reference end are excluded (`:527-530`).
- Records are written only if ≥1 true SV exists, or the signal was `LRS_SIG_X` (`SV_core_SRS_ins_del.cpp:3992`).
- Labels: net length difference <10 is "INDEL"; otherwise "DEL" or "INS"; single bases are "SNP" (`:4025-4071`).
- **No merging of adjacent gaps:** the summed-INS/DEL/CPX code is disabled with `if(false)` (`Contig_aligner.hpp:415, 434, 454`; that is the SRS function, and the TGS function has no summing). One SV can therefore split into several records.
- Filters when support ≤1:
  - >5 SVs ≥ MIN_sv_len: drop all (`:4094-4095`);
  - 2 to 5 SVs: opposite-sign SVs <10 bp apart drop all (ERR_PRONE: any two <100 bp apart) (`:4098-4136`);
  - support ==1 and an ALT whose 10-mers recur ≥6 times make up >40%, with no shared 10-mer in the ±100 bp flanks: drop all (`:4142-4241`).

## 6. Contig polishing with SRS (`SV_core_SRS_ins_del.cpp:3191-3213`; `Contig_polishing.hpp`)

- **When it runs:** an SRS BAM is present and a `-k` profile is loaded (`kmer_error_profile != NULL`) (`:3192`), and the contig is longer than 50 bp (`:3208`). It applies to every selected haplotype.
- **Profile** (`RefHandler.hpp:340-365`): a TSV of k-mer and probability, with a header line skipped. k comes from the first line's k-mer length. The table has 4^k floats, so a 9-mer profile holds 262,144 entries (`:353-355`). The code works for any k; 9 is not hardcoded.
- **Which positions get polished** (`RefHandler.hpp:370-402`):
  - a contig position is 'P' if its profile k-mer is all-A or all-T, or has profile value <0.6;
  - positions within 25 k-mers of a 'P' become 'E', i.e. indexed as anchors (`:364, 392-401`).
- **Algorithm** (`Contig_polishing.hpp:508-517`):
  1. Index contig 13-mers at non-'N' positions (`:195, 202-224`).
  2. Fetch SRS reads in [region ±300] that are not duplicate, secondary or supplementary, with MAPQ ≥20 and length <1024 (xassert) (`:267-283`; `SV_core_SRS_ins_del.cpp:3199-3202`).
  3. Anchor each read by its best diagonal with >25 unique 13-mer hits (`Contig_polishing.hpp:226-265`).
  4. Collect each read's sequence over 'P' runs (`:319-329`). Where one string has ≥3 reads and ≥70% of reads, replace the contig segment if it differs by 1-2 bases (`:116-127, 476-495`).
  5. Homopolymers ≥8 in reads are clustered within 10 bp. If one length has >70% and ≥3 reads, resize the contig run (`:332-344, 378-474`).
- **Pitfall:** polishing edits the member `bin_contig` only (`SV_core_SRS_ins_del.cpp:3206-3210`). `LRS_SV_generating_Germline` rebuilds `bin_contig` from the unpolished `strContig` (`:3366`). Somatic rebuilds `bin_tc` (`:3650`). **The polished sequence is never used downstream** (no write-back is present). See pitfalls-B.
- `use_contig_num` limits the loop to final paths (`:3193`).

## 7. Genotyping

### LRS (`SV_core_LRS.cpp:2279-2312`; `SV_core.hpp:1106-1113`; `NovaSVRst.hpp:595-597, 622-665`)
- **Support matrix:** `supp` is the haplotype's full-cover or clip best-target read count. `not_supp` is the other final haplotype's support when 2 haplotypes exist, else 0 (`SV_core.hpp:1107-1109`). There is no per-read REF realignment; reads assigned to no final contig are not counted as REF.
- **GT** (`calGT_LRS_CORE`, `SV_core_LRS.cpp:2279-2290`):
  | Condition | GT |
  |---|---|
  | not_supp > 6·supp | 0/0 |
  | 6·not_supp < supp | 1/1 (0/1 when supp==1) |
  | otherwise | 0/1 |

  Two 0/1 records on different haplotypes of one region become 1/1 (`SV_core_LRS.cpp:2444-2465`; `NovaSVRst.hpp:1833-1844`).
- **QUAL/GQ** (`get_SV_quality_score_core`, `SV_core_LRS.cpp:2292-2312`). Pseudo log-likelihoods:
  - RR = −5·alt − 5.3
  - AA = −5·ref
  - RA = −0.3·(ref+alt) − 5
  - QUAL = 10·(max(RA,AA) − RR), clamped to [20, 500] (`NovaSVRst.hpp:565-566`). It is divided by 10 if unknown/alt > 0.1 and QUAL > 100.
  - In VCF output, `unknown_read_n` is passed as 0 (`SV_core_SRS_ins_del.cpp:3415`).
- **FILTER** (`NovaSVRst.hpp:660-665`):
  - |SVLEN| < 50: always PASS;
  - else 0/0: LOW_DEPTH. **0/0 records are never written** (`writeVCF_final` returns false, `:624`);
  - for LRS calls the remaining checks (`DR_not_support`, VNTR LOW_QUAL) default to pass (`NovaSVRst.hpp:385`).
  - In practice every emitted LRS record is PASS, and QUAL is the only quality signal.
- `GT_cal.hpp`'s Bayesian `GT_handler` (priors at `GT_cal.hpp:130-160`) is **dead code**: no reference exists outside the file, and `GT_cal.cpp` is empty (`GT_cal.cpp:1-10`).

### SRS re-GT of LRS calls (hybrid)
- `getGT_from_support_read` (`NovaSVRst.hpp:1507-1528`):
  - alt > 4·ref: 1/1;
  - 3·alt < ref: 0/0;
  - else 0/1;
  - unknown > 3·(alt+ref), or alt < 2: 0/0.
- Applied under the conditions in control-flow §1c (`SV_core.hpp:2523-2549`).

## 8. Somatic mode (`-X` tumor, `-x` normal; "only supports LRS (ONT or HiFi)" `sv_call_main.hpp:111-113`)

1. Tumor signals: hard clips >1000 bp become clip signals too (`SV_core_LRS.cpp:2119-2132`). The support gate is `depth_T*0.1` (`SV_core.hpp:2371, 2394`).
2. Tumor and normal are assembled separately with `loading_SUPP_READ=true`: reads ≥100 bp, keeping double-clipped reads (`SV_core_SRS_ins_del.cpp:1009-1011`). The target cap is 15 (`:3060-3062`).
3. Per tumor contig (`SV_core_SRS_ins_del.cpp:3591-3870`), with `MIN_READ_N_TUMOR=3` and `MIN_READ_N_NORMAL=1` (`:3599-3600`):
   - **Full-cover** (≥3 full-cover reads, or ≥3 left- and ≥3 right-clip reads) (`:3629-3632`):
     - compare with each full-cover normal contig: length difference must be <max(15%, 50 bp) (`:3461-3466`);
     - then ksw2 T-vs-N with the somatic preset, z-drop = max len, band = 0.3·max len (`:3481-3484`);
     - "shared" if the alignment is end-to-end with no indel > min_SV_diff (50; 200 or 500 for long low-support normal contigs) (`:3472-3523`);
     - shared contigs are skipped; others are realigned to the reference like germline calls and tagged TUMOR_UNIQ, or UNKNOWN if the normal list is empty (`:3655-3716`).
   - **Long-range** (≥3 clip reads) (`:3633-3634, 3719-3855`):
     - shared if the normal has a same-side clip contig;
     - otherwise SA tags of supporting reads give BND breakpoint pairs, clustered within 20 bp on both ends (`SV_core.hpp:1380-1449`);
     - the largest cluster with ≥2 reads emits a BND, plus a reciprocal BND if same chromosome and <1000 bp (`SV_core_SRS_ins_del.cpp:3801-3830`).
4. GT uses the same `calGT_LRS`. BND records get `not_supp=0` (`:3580-3581`). No SRS re-GT is done.

## 9. Limits relevant to non-human genomes and satellite arrays

| Limit | Value | Line |
|---|---|---|
| Block size | 5 Mb, non-overlapping | `RefHandler.hpp:30, 248` |
| Whole chromosome held in RAM (uint8 per base) | chromosome length in bytes | `RefHandler.hpp:314-321` |
| Max reads per assembly region | `MIN(MAX(5·depth,100),1000)` germline; tumor `MAX(10·depth,100)` | `SV_core_LRS.cpp:2320-2321`; `SV_core.hpp:2378` |
| Read soft-clip slice | clips ≥100 kb are not loaded | `ReadHandler.hpp:810` |
| Long-range SV size from ref-index clustering | 1 kb < SV < 50 kb | `SV_core.hpp:1606-1607` |
| Ref index uniqueness | k=16 sampled every 3 bp, unique only; **repeats/satellites produce no anchors** | `RefHandler.hpp:209, 283, 297` |
| Contig k-mer anchor | unitigs occurring >30× in a target path ignored | `SV_core_LRS.cpp:630` |
| Long alignment padding | ≤20 kb (10 kb per side); band/zdrop 30 kb | `Contig_aligner.hpp:181-205` |
| Single-support contigs >10 kb in clip/long-range mode | skipped | `SV_core_SRS_ins_del.cpp:3392-3395` |
| Contig length in PURE_STR | ≤50 kb (3 Mb in ASM) | `SV_core_LRS.cpp:2324-2328` |
| VCF line buffer | 1,000,000 chars, fixed (`sprintf`, no resize) | `SV_core.hpp:2312-2313`; `NovaSVRst.hpp:628-669` |
| `<DEL>` symbolic | DEL < −50 kb | `NovaSVRst.hpp:643` |
| Ploidy | diploid hardcoded: ≤2 final contigs, 0/0, 0/1, 1/1 only, homozygous merge | `SV_core_LRS.cpp:1265-1269, 2279-2290, 2444-2465` |
| SA supplementary min length | 1200 bp | `SV_core_LRS.cpp:1999` |
| ASM preset | 1 sequence per region (haploid), stitched ≤3 Mb | `SV_core_SRS_ins_del.cpp:1029, 1253-1274` |
