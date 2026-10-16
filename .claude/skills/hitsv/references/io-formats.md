# Input and output formats

Scope: HitSV v2.0.4, commit 93880367. Paths are relative to `src/` unless they start with `demo/`. All line numbers were verified by reading the source.

## 1. Inputs for `call`

Options are registered in `SVcalling_core/sv_call_main.hpp:96-245`.

| Flag | Input | Required | Cite |
|---|---|---|---|
| `-r` | reference FASTA; `.fai` is built if missing | always | 105, 188-191; `clib/bam_file.c:262-274` |
| `-l` / `-n` | LRS BAM / SRS BAM (indexed: .bai/.csi/.crai) | at least one of `-l`, `-n`, `-X` | 107-108, 192-195; `clib/bam_file.c:81-93` |
| `-X` / `-x` | tumour / normal LRS BAM (somatic) | optional | 112-113 |
| `-I` | local-repetitiveness file from `srs_fa_stat` | **mandatory when `-n` is given** | 116, 197-200 |
| `-T` | SRS BAM-stat JSON from `ngs_bam_stat` | optional; computed on the fly if absent | 158, 273-277 |
| `-L` | SRS translocation-read BAM from `srs_trans_reads` | optional | 157 |
| `-R` | realigned SRS reads (needs an `SV:Z:` tag) | optional; loader disabled by `if(false)` | 159; `SV_core.hpp:2590`; `SV_core_SRS_ins_del.cpp:524-526` |
| `-k` | LRS k-mer error profile (`kmer\tprob`) | optional (warning for hybrid) | 163, 259-262, 280-282 |
| `-b` | BED of forced-contig regions (LRS only) | optional | 169, 287-289 |

### 1.1 Read filters

#### SRS block scan (`SRS_GET_SIGNALS_FROM_READs`, `SV_core_SRS_ins_del.cpp:419-501`)

| Step | Rule | Cite |
|---|---|---|
| 1 | drop duplicates (before depth counting) | 437-438 |
| 2 | drop secondary and supplementary alignments | 444-447 |
| 3 | an unmapped read becomes a UM signal (`pass_compact_filter`: ≥ 8 runs after homopolymer compression) | 449-452, 365-375; `ReadHandler.cpp:11-20` |
| 4 | MAPQ: pass if `MAPQ ≥ 20`, **or** mate `MQ` tag ≥ 20 with the same chromosome and \|isize\| < 20000 | `SV_core_SRS_ins_del.cpp:17-25`; MQ read from tag at `clib/bam_file.c:565` |
| 5 | split/clip read (SR) if a clip exists or `max(gaps, NM) ≥ 2` | 393-417 |
| 6 | SR storage filters: compact ≥ 8; reject if ≥ 19 of the 20 bases next to a clip are A or T (the code tests `> 18`); clip average quality < 0.6 × mapped average quality drops that clip | `ReadHandler.cpp:78-92`, 23-45, 59-76 |
| 7 | DR (discordant pair): mate mapped and (different chromosome, same strand, isize ≤ min or ≥ max, or wrong orientation) | `clib/bam_file.c:509-530` |
| 8 | translocation-read list: different chromosome or \|pos − mpos\| > 200000 | 471-473 |

Per-block depth estimate: `read_counter * MaxReadLen / block_len` (`:476`). `read_counter` counts every non-duplicate read, including secondary and unmapped ones (`:437-442`).

The SRS assembly reload (`assembly_load_read`) drops duplicate, secondary and supplementary reads (`SV_core_SRS_ins_del.cpp:1325-1327`).

#### LRS block scan (`LRS_SIG_COLLECTION`, `SV_core_LRS.cpp:1945-2178`)

| Rule | Cite |
|---|---|
| drop secondary | 1994-1995 |
| drop MAPQ < 1 (MAPQ 0) | 1996-1997 |
| drop supplementary reads with `l_qseq < 1200` | 1999-2000 |
| **duplicates are not filtered** (no `bam_is_duplicate` in this file) | grep → 0 hits |
| error-rate check (skipped for `-p ASM`): per 1000 aligned bases, mismatch + indel-event rate > 0.2 rejects the read | 2004, 2468-2512 |
| whole-read rate > **0.04** rejects the read (> **0.16** with `ERR_PRONE`); hard-coded, `-e` is ignored | 2517-2526 |
| at assembly: drop supplementary (0x800) reads with MAPQ < 20; drop reads clipped > 200 bp on both sides; drop reads shorter than `min_LRS_read_length` (somatic: < 100 bp) | `SV_core_SRS_ins_del.cpp:1005-1017` |

LRS depth estimate: `Σ primary_len / (mean_len + block_len)` (`SV_core_LRS.cpp:2166-2173`).

## 2. Reference blocks and region options

- **Block size:** `SEGMENT_LEN 5000000` (`RefHandler.hpp:30`). The loop comment says "2M bases per block", which is stale (`sv_call_main.hpp:298`).
- **Algorithm** (`RefHandler.hpp:243-267`):
  - `init` sets `cur_r = [st_pos - 5M, st_pos]` on `st_chr_ID` (`:202`).
  - Each call advances `st` by 5 M and sets `ed = st + 5M - 1` (`:248`).
  - If the new start (`old_st + 5M`) would pass `len-1`, it jumps to `chr+1, [0, 5M]` (`:245-246`). That block's end is 5,000,000 inclusive, so it is 1 bp longer than the others. Every block end is then clamped to `len-1` (`:250`).
  - Consecutive blocks do **not** overlap, except for 1 bp between the first and second block of each chromosome reached by the jump: [0, 5M] and then [5M, 10M−1]. Example (traced by hand): a 12 Mb chromosome gives [0,5M], [5M,10M−1], [10M,12M−1], and then the jump. There is no flanking overlap. Calls near block edges rely on the assembly load margins (`region_addition_load` 500–50000, `SV_core.hpp:2230-2233`).
  - Calls that overlap a block boundary are kept if their interval overlaps the block (`SV_core.hpp:2389-2391`, `SV_core_LRS.cpp:2334-2336`).
- **`-S/-E`** are 0-based **BAM header indices**, not names (`sv_call_main.hpp:125-126`; defaults 0 and 50,000,000).
- **`-s/-F`** (defaults 0 and 1e9, `:127-128`): `-s` applies only to `-S`, through the initial `cur_r` (`RefHandler.hpp:202`). `-F` applies only when `chr == ed_chr_ID` (`:251-252`). Iteration stops when `chr > ed_chr_ID`, or when `chr == ed_chr_ID` and `st > ed_pos` (`:255-258`).
- With `-s 100000` the first block is `[100000, 5099999]`. Blocks align to `-s`, not to multiples of 5 Mb.
- The reference is fetched by the BAM header contig name (`RefHandler.hpp:314`), so BAM and FASTA names must match. The whole chromosome is held in memory as a DNA5 code per base (`:315-321`).

## 3. Local repetitiveness file (`srs_fa_stat` → `-I`)

### 3.1 Algorithm (`SVcalling_core/fa_stat.hpp:28-106`)

1. Run as `HitSV srs_fa_stat ref.fa > ref.stat.txt` (`main.cpp:88-91, 108-109`). `run_ana` requires exactly one argument (`fa_stat.hpp:112-119`).
2. Parameters: window `REGION_LEN = 600`, step `REGION_STEP_LEN = 300`, k from `MIN_K_size = 20` to `MAX_K_size = 180` in steps of `K_step_LEN = 20` (`:34-38`). Every FASTA record is processed, with no chromosome cap (the `MIN(26, N_seq)` cap is commented out, `:44-45`).
3. For each window `[begin, min(begin+600, chrLen))`, with `begin += 300` (`:59-61`), skip windows shorter than 180 bp (`:62-63`).
4. For k = 20, 40, …, 180: insert every k-mer that contains no `'N'` into a `std::set<std::string>`. On the **first repeated k-mer**, set `withloop = true` and `kmerSizeF = k`, then try the next k (`:70-87`). The first k that has **no** repeat ends the loop (`:89-91`).
5. Result R = **the largest k in {20…180} for which some k-mer occurs at least twice in the window**. That is the longest repeated k-mer, rounded down to a multiple of 20 and capped at 180. A window with no repeated 20-mer prints nothing (`:95-101`).
   - Caveat: once one k has repeated, `withloop` is never reset to false (`:68`, only set at `:84`). So for k > 20 the loop cannot break and runs on to k = 180, and `kmerSizeF` ends at the **last** k that had a repeat. This equals the longest repeated k for real sequence, because a repeated k-mer implies a repeated shorter k-mer. Cost: the loop always runs 9 passes once a 20-mer repeats.
   - `N_count` checks only uppercase `'N'` (`:17-24`). Lowercase (soft-masked) bases are compared case-sensitively (htslib keeps case; `htslib/faidx.c:712` copies all printable characters). So `acgt` and `ACGT` count as different k-mers.

### 3.2 File format

Tab-separated with no comments:

| Line | Format | Cite |
|---|---|---|
| header (1st line) | `600\t300\t0` (window, step, unused) | `fa_stat.hpp:48-50` |
| data | `<chrIdx>\t<windowStart>\t<R>`; chrIdx is the 0-based **FASTA record index** (`faidx_iseq(i)`), windowStart is 0-based | `:97` |
| BED variant (`run_ana_bed`, `:123-134`) | `chrName\tstart\tend\tR` | `:99` |

The BED variant is not reachable from `main` (`main.cpp:105-115`).

Example from `demo/GRCh38.stat.txt.gz`: the header is `600 300 0` and the first rows are `0 9600 100`, `0 9900 100`, `0 10200 140`. It covers 25 chromosome indices. The R histogram is 20: 1,048,217; 40: 141,898; 60: 28,977; 80: 10,980; 100: 6,132; 120: 3,558; 140: 2,272; 160: 1,461; 180: 4,051.

### 3.3 Loading and querying (`RefHandler.hpp:102-141`)

- Every line must have exactly 3 tab-separated fields, or `xassert` aborts (`:111`). The first line sets `region_length` and `region_step_size` (`:112-116`).
- Only rows with **R ≥ 40** are stored (`:123-125`), so the 20-mer-only rows (about 84% of the GRCh38 file) are ignored. Key = `pos + (chrIdx << 32)` (`:120`).
- `search(chr, bg, ed)` probes each integer position from `bg - bg%300` to `ed - ed%300 + 300` (`:132-139`). Only window starts that are multiples of the step can hit, so a file built with another step would mostly miss.
- Uses:
  - Every hit in the current 5 Mb block becomes a forced VNTR INS signal `SVE(chr, start+300, ..., score 99, SOLID, range 300)` (`SV_core_SRS_ins_del.cpp:663-671`, `get_middle` at `RefHandler.hpp:93-96`). For GRCh38 that is 199,329 rows with R ≥ 40, i.e. about 199k forced VNTR candidates genome-wide (sum of the histogram in §3.2). Each one is later reduced by `SRS_combine_duplication` (`:675-676`).
  - Small-variant/VNTR step: `rlr_result` for each signal is printed in the log. The VNTR caller is triggered by `withRepeat` from assembly, not by `rlr` (`SV_core_SRS.cpp:537-556`).
  - Hybrid: LRS calls in a repeat window are not re-genotyped when depth < 10 and `-p ERR_PRONE` (`SV_core.hpp:2529-2536`).
- **Chromosome identity:** chrIdx is the FASTA order, while the query uses the BAM `tid` (`SV_core_SRS_ins_del.cpp:664`). The file is only valid if the FASTA record order equals the BAM header order.
- **Missing chromosomes:** a chromosome with no rows simply gets no hits (`find` misses, `:136-139`). No error, and SRS calling still runs, without VNTR seeding there.
- **Missing file:** the `else` path that tolerates an absent `-I` (`sv_call_main.hpp:267-271`) cannot be reached, because `get_option` already rejected `-n` without `-I` (`:197-200`).
- **Hybrid crash path:** `xassert(rlr->region_length != 0)` in the hybrid re-genotyping step (`SV_core.hpp:2530`) aborts if `-I` was an empty or headerless file.

**Building it for a non-human genome:** run `HitSV srs_fa_stat genome.fa > genome.stat.txt` on the exact FASTA used for alignment (same record order as the BAM header). Nothing in the algorithm is human-specific. Runtime is roughly proportional to the number of windows times up to 9 `std::set<string>` passes of about 580 strings each (inferred).

## 4. SRS BAM-stat JSON (`ngs_bam_stat` → `-T`)

Run as `HitSV ngs_bam_stat ref.fa srs.bam > stat.json` (`main.cpp:93-98, 106-107`). The help text says `srs_bam_stat`, but the registered command is `ngs_bam_stat`.

| Key | Meaning | How computed | Cite |
|---|---|---|---|
| `minInsertLen` | 1st percentile isize | `getInsertLen(…, 0.01)` | `bam_stat.hpp:126` |
| `middleInsertLen` | median isize | `0.5` | 127 |
| `maxInsertLen` | 99th percentile isize | `0.99` | 128 |
| `isize_distribution` | array, one entry per isize in `[min, max)`: count / (high-confidence pairs + 1) | 130-137 |
| `analysis_read_length` | the first read length bin holding > 60% of reads; falls back to the mean | 105-118 |
| `ave_read_len` | mean read length over the first 100,000 primary non-dup reads (lengths < 1000 only enter the sum) | 72-116 |
| `ave_read_depth` | mean of per-slice depth = bases / span over sampled slices | 125; `cpp_lib/statistics/StatsManager.cpp:136-205` |

Further details:
- Read-length pass: skips secondary, supplementary and duplicate reads, and stops after 100k reads (`bam_stat.hpp:88-96`).
- isize pass: reads are filtered by `CoreReadFilter`, which requires a mapped pair, MAPQ > 0, no SA split, and a simple alignment; one pair per start position (`StatsManager.cpp:13-92`). Sampling stops when the insert-size distribution converges or after 5 M observations (`StatsManager.hpp:109-117`).
- `json_load` reads the same keys (`bam_stat.hpp:26-49`). Lines 42-44 re-read `minInsertLen` three times (harmless).
- **Without `-T`**, `call` runs `generate_state(ref, srs_bam, false)` (`sv_call_main.hpp:276-277`). This does the same two passes over the whole BAM, at the start of every run including each `-S/-E` region job, and prints nothing (`bam_stat.hpp:140-157`).

How the stats are used downstream:
- `insert_size_min/max` (`SV_core.hpp:61-65`);
- the DR/SH/UM breakpoint distributions (`bam_stat.hpp:159-217`);
- the depth-counter scale (`SV_core.hpp:2267`).

`getBreakPoint_Distribution` indexes `isize_distribution[i + 2*readlen - min]` (`bam_stat.hpp:173`). It needs `maxInsertLen > 2*analysis_read_length`, otherwise the vectors are empty or negative-sized (inferred).

## 5. VCF output

### 5.1 Header (`NovaSVRst.hpp:1638-1697`)

The header is written in `init` unless `-H` or PURE_STR is set (`SV_core.hpp:2314-2315`).

| Kind | IDs | Cite |
|---|---|---|
| meta | `##fileformat=VCFv4.1`, `##fileDate`, `##source=HitSV2.0.0` (`PACKAGE_VERSION "2.0.0"` in `clib/desc.h:2`) | 1644-1646 |
| INFO | SVTYPE, BND_TYPE, SVLEN, END, IMPRECISE, VNTR, VNTR_DepthLow, CONTIG, CIGAR, SIGNAL_REGION, CONTIG_POSITION, SOMATIC | 1649-1660 |
| FORMAT | GT, GQ, SR, LR, SL, IF | 1663-1670 |
| FILTER | LOW_DEPTH, lt50bp, DR_not_support, Duplication_sv, LOW_QUAL | 1672-1676 |
| ALT | DEL, INS, SNP, INDEL, TRA, INV, BND | 1679-1685 |
| contig | one `##contig=<ID,length>` per BAM target | 1688-1689 |
| command | `##HitSV_call_Command=<argv>; Date=<ctime>` and `##HitSV Version=…` | 1692-1694; `sv_call_main.hpp:64-94` |
| columns | `#CHROM … FORMAT Sample` (sample column always named "Sample") | 1696 |

The `##fileDate` value is computed as `1990 + tm_year`, with month and day printed without zero padding (`:1645`). This explains the demo values `2115116` and `211529`.

### 5.2 Record (`writeVCF_final`, `NovaSVRst.hpp:622-746`)

- **Skipped entirely** when GT is 0/0 or 0|0, or `is_duplication_sv` (`:624`). So `LOW_DEPTH` and `Duplication_sv` are effectively never printed.
- CHROM is the BAM target name. POS is `max(st_pos, 1)` (`:628`).
- ID = `<globalID>_<chr>_<SVTYPE>_<st>_<end|0>_<|len|>` (`:639-640`).
- REF/ALT: full sequences, except non-BND deletions longer than 50 kb, which are written as `<first base>\t<DEL>` (`:643-644`).
- QUAL: see 5.3.
- FILTER, in order (`:660-665`):
  1. \|SVLEN\| < 50 → `PASS`. The header-declared `lt50bp` is never produced by this code, although `demo/lrs.vcf` contains it (written by an older build, `##source=GC_SV2.1.0`).
  2. otherwise `Duplication_sv`, `LOW_DEPTH` (GT==0), `DR_not_support` (`!deletion_pass_DR_filter`), `LOW_QUAL` (VNTR caller with `QUAL_presetVNTR ≤ 10`), or else `PASS`.
- INFO:
  - BND records: `[IMPRECISE;]SVTYPE=BND;BND_TYPE=…` (`:668-671`).
  - Other records: `SVTYPE=…;END=…;SVLEN=…[;IMPRECISE][;VNTR[;VNTR_DepthLow]]` (`:672-684`).
  - Always `;SIGNAL_REGION=chr:st-ed` (`:687`).
  - If there is a CIGAR and \|SVLEN\| ≥ 10: `;CONTIG_POSITION=chr:pos;CIGAR=…;CONTIG=<seq>` (`:690-703`).
  - `;SOMATIC=…` for somatic calls (`:705-708`).
- FORMAT is always `GT:GQ:LR:SR:SL:IF` (`:713`):
  - `LR = not_supp,supp,unknown` (long reads); `SR = REF,ALT,unknown` (short reads, `[2][1],[2][0],[2][2]`) (`:736-738`).
  - Pure SRS calls get LR = `-1,-1,-1` and IF = `-1×5` (`:727-728`). LRS calls with no SRS data get SR = `-1,-1,-1` (`:729-733`).
  - `SL` = assembler-suggested length (`:740`). `IF = regionID,hapID,totalHap,SVinHap,totalSVinHap` (`:742`).

**SVTYPE values emitted:**
- `INS` and `DEL` (`Contig_aligner.hpp:373,399,429,448,474,487`; `SV_core_SRS_ins_del.cpp:4031-4050`).
- `INDEL` when the size difference is < 10 (`SV_core_SRS_ins_del.cpp:4025-4029, 4043-4046`) and `SNP` (`:4068`); these are LRS small variants.
- `INV` from SRS (`SV_core_SRS.cpp:402`).
- `BND`: SRS TRA via `store_BND` with `BND_TYPE` = `SV::STR`, e.g. TRA (`SV_core_SRS_inv_tra.cpp:202`); somatic LRS BND with `BND_TYPE=NA` (`SV_core_SRS_ins_del.cpp:3583-3587`).
- **No DUP and no TRA SVTYPE** is ever emitted (TRA appears only as BND_TYPE). Only precise TRA BNDs are written (`SV_core_SRS.cpp:455-470`).

### 5.3 Genotype and quality

- **LRS GT** (`calGT_LRS_CORE`, `SV_core_LRS.cpp:2279-2290`), with R = not_supp and A = supp:
  - 0/0 if R > 6A;
  - 1/1 if 6R < A, except A == 1 gives 0/1;
  - 0/1 otherwise.
- Haplotype merge: identical alleles from different haplotypes of the same region are merged (supp summed, not_supp reduced), which turns 0/1 + 0/1 into 1/1 (`SV_core_LRS.cpp:2444-2458`, `NovaSVRst.hpp:1840-1844`).
- **SRS GT** (`getGT_from_support_read`, `NovaSVRst.hpp:1507-1528`), the same rule for INS and DEL with A = ALT reads and R = REF reads:
  - 1/1 if A > 4R;
  - 0/0 if 3A < R;
  - else 0/1;
  - forced to 0/0 if 3(A+R) < unknown or A < 2.
- Extra SRS deletion rules (`NovaSVRst.hpp:1592-1635`):
  - DR filter for VNTR deletions < −100 bp and for deletions < −400 bp;
  - DR rescue for 0/0 deletions < −700 bp;
  - deletions > 50 kb are forced to 0/0, so they are dropped;
  - read-depth filter for deletions > 1 kb (ratio > 0.9 → 0/0).
- Force-ALU calls halve the REF count (`:1577-1578`).
- Hybrid: an LRS call is set to 0/0 when supp == 1, not_supp ≥ 1 and SRS REF ≥ 20 × ALT. It is re-genotyped by SRS when not_supp == 0 and supp < 4 (`SV_core.hpp:2545-2548`, `NovaSVRst.hpp:599-613`).
- **QUAL/GQ** (`get_SV_quality_score_core`, `SV_core_LRS.cpp:2292-2312`), with R = REF reads, A = ALT reads, U = both/unknown:
  - log-likelihoods: `RR = −5A − 5.3`, `AA = −5R`, `RA = −0.3(R+A) − 5`;
  - `QUAL = 10·(max(RA, AA) − RR)`, clamped to [20, 500]; divided by 10 if U/A > 0.1 and QUAL > 100; set to 20 if the length is imprecise;
  - `GQ = 10·(L_GT − max(other two))`, clamped to [20, 500];
  - VNTR caller results use `QUAL = GQ = QUAL_presetVNTR` (`NovaSVRst.hpp:653-655`).
- **GT_cal.hpp is dead code.** `class GT_handler` (`GT_cal.hpp:42`) is only used by `GT_CAL_DEMO_10_10_10_10` inside the same header (`:243-247`). GT_cal.cpp is empty (`GT_cal.cpp:1-10`). It is only `#include`d by `main.cpp:11`.

### 5.4 Random phasing (`-f`, default **1 = on**, `sv_call_main.hpp:133`)

`randomly_phasing_Genotype` (`NovaSVRst.hpp:520-537`):
- 0/0 → 0|0, and 1/1 → 1|1.
- 0/1 → 0|1 or 1|0:
  - if the region had more than one haplotype, hap 0 gives 0|1 and the others give 1|0 (local phasing);
  - otherwise the choice is `random_num % 2`.

The "random" number is `(int)(global_SV_ID * 1.897239)` (`SV_core.hpp:2453,2485,2556,2619`), which is deterministic. `global_SV_ID` only increases inside `writeVCF_final` (`NovaSVRst.hpp:632-634`). In the somatic branch, the basic LRS branch and the SRS branch, all calls are phased in one loop before any is written (`SV_core.hpp:2451-2459, 2554-2564, 2617-2631`). So every single-haplotype 0/1 in a block gets the **same** orientation; only the FC-BED branch interleaves phasing and writing (`:2483-2488`). SRS calls are phased with the same expression (`SV_core.hpp:2617-2621`), so SRS-only VCFs show `|` genotypes by default (demo/srs.vcf: `1|0`). An unknown GT code would hit `xassert(0)` (`:535`).

### 5.5 Small-variant output (`-v`, default 1)

- It is forced off for `-p ERR_PRONE` (`sv_call_main.hpp:135, 219-220`).
- With `-v` on, `get_canditate_SVs_TGS` records every mismatch inside M runs and every I/D (excluding the first and last CIGAR op, and insertions within 20 bp of the region end) (`Contig_aligner.hpp:508-553`).
- These become SNP / INDEL / INS / DEL records with LRS_INFO (`SV_core_SRS_ins_del.cpp:4019-4073`), alongside the SVs of the same contig. This only happens if the contig carries at least one variant of `MIN_sv_len` or more (`truSV_num != 0`) or is a replacement (`:3990-3992`).
- They go into the same VCF, so default LRS output mixes SNPs and indels with SVs. The `-m` value (default 25) is only a calling threshold, not an output filter for these records.

## 6. PURE_STR output (`-B PURE_STR`)

- No VCF header (`SV_core.hpp:2314`) and no k-mer index build (`sv_call_main.hpp:304-305`).
- **Only produced in the `-b` FC_BED branch.** The basic germline path calls `LRS_SV_CALLING_germline(..., false)` (`SV_core.hpp:2506`), so `-B PURE_STR` without `-b` writes nothing.
- One line per signal region to the output (`SV_core_LRS.cpp:2338-2407`):
  - `REF_REGION=chr:POS-END\tSTAT=NO_ENOUGH_READ_SUPPORT` (`:2340-2341`);
  - `…\tSTAT=ASSEMBLY_FAIL` (`:2349-2350`);
  - `…\tSTAT=NO_CONTIG_SELECTED` (`:2376-2377`);
  - for each PASS haplotype contig: `REF_REGION=chr:loadSt-loadEd\tSTAT=PASS\tGT:QU:GQ:GI:CI:CN=<GT>:<QUAL>:<GQ>:<GTint>:<contigIdx>:<nContigs>\tSUP_R:UNSUP_R:FULL_COVER_R:CLIP_AT_LEFT:CLIP_AT_RIGHT:IS_FINAL=…\tCONTIG_STR=<seq>` (`:2398-2406`).
- QUAL is computed with "imprecise" = true, so QU is always 20 (`:2395` → `SV_core_LRS.cpp:2303-2304`).
- Contigs over 50 kb are skipped, or over 3 Mb with `-p ASM` (`:2324-2328, 2367`).
- Any `-B` value other than `VCF`/`PURE_STR` leaves `OUT_MODE_UNSET`. VCF records are then still written in the basic path (it does not test `OUT_MODE`), but the k-mer index is not built (`sv_call_main.hpp:237-242, 304`).

## 7. Other outputs

- **stderr** gets very large logs: every SRS VCF line is echoed (`SV_core.hpp:2625-2630`), and so is every LRS line (`:2560-2563`). The README recommends `2> /dev/null`.
- `##contig` lines and CHROM use the BAM header of the last opened BAM (`SV_core.hpp:2278-2292`). SRS INV/BND are written with the SRS header (`SV_core_SRS.cpp:461, 611`).
