---
title: "hitsv: Architecture"
type: reference
generated: true
---

> Part of the [hitsv](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/hitsv/references/architecture.md)

{% raw %}
# Architecture

Full commit SHA: `93880367ef76f6440c918f7dd9744dda2659bcdb` (2026-08-12, "Add files via upload";
the upstream repo is a single squashed commit, so there is no history to bisect).

Paths are relative to the repo root.

## Directory map

| Path | Role | Key files |
|---|---|---|
| `src/main.cpp` | Subcommand dispatcher (118 lines) | `main()` `:100-117` |
| `src/SVcalling_core/` | Everything in `call`: options, reference blocks, read loading, signals, assembly driver, VCF writer | `sv_call_main.hpp`, `SV_core.hpp`, `SV_core_LRS.cpp`, `SV_core_SRS*.cpp`, `NovaSVRst.hpp` |
| `src/SVcalling_core/ANALYSIS/` | `HitSV tools` sub-tools: cohort CSV windows, contig clustering, VNTR / genome-context stats | `analysis.cpp:160-185`, `FC_joint_calling.hpp` |
| `src/cpp_lib/Assembler/` | Local assemblers | `DBGAssembler.hpp` (`MainDBGHandler1` `:3095`), `GreedyAssembler.*` (SRS only) |
| `src/cpp_lib/` | Option parser, JSON (cJSON), stats, graph utilities | `get_option_cpp.hpp` (`COMMAND_HANDLER` `:332-384`) |
| `src/clib/` | C helpers for BAM/VCF records | `bam_file.c` (DRP orientation/insert tests), `utils.c` |
| `src/kswlib/` | ksw2 `extd2_sse` used for contig↔reference realignment | `ksw2_extd2_sse.c` |
| `src/htslib/` | Vendored htslib (statically linked) | — |
| `py/` | Python post-processing for cohort CSV / MEI-TR arrays | `mei_tr_array_single_region_detection.py` |
| `HitSV_multy_process.sh` | Per-chromosome parallel wrapper (xargs -P) | see pitfalls #1 |
| `cohort_csv_analysis.md` | Author's cohort CSV recipe (contains broken commands, see `pitfalls.md` D) | — |
| `demo/` | 1 Mb hs37d5 slice, 10× HiFi BAM, `9mer_profile`, GRCh38 stat file | README calls it `demo_HitSV/`; SRS demo BAM and `SRS_HG002_stat.json` are missing |
| `Release/` | Eclipse-generated makefile **plus 70 tracked `.o/.d` files** and two prebuilt binaries | run `make clean` first |

## Entry points

| Command | Entry | Notes |
|---|---|---|
| `HitSV call` | `main.cpp:36-44` → `NOVA_SV_PARA::get_option` `sv_call_main.hpp:96-245` → `NOVA_SV_handler::run` `:252-312` | Loops 5 Mb blocks (`RefHandler.hpp:30`), calls `SV_CALLING_Handler::SV_calling_in_region_main_pipeline()` `SV_core.hpp:2355` |
| `HitSV ngs_bam_stat ref.fa srs.bam > stat.json` | `main.cpp:93-98` → `bam_stat.hpp` | Help text wrongly says `srs_bam_stat` (`main.cpp:107`) |
| `HitSV srs_fa_stat ref.fa > ref.stat.txt` | `main.cpp:88-91` → `fa_stat.hpp:34-101` | Local-repetitiveness R file for `-I` |
| `HitSV srs_trans_reads ref.fa in.bam out.bam` | `main.cpp:46-86` | Keeps primary, non-dup pairs with \|isize\| > 5000 or 0; swaps pos/mpos so mates can be fetched by region |
| `HitSV ngs_merging`, `HitSV ngs_fc` | `main.cpp:24-34` | Developer leftovers with hard-coded paths; unusable |
| `HitSV tools <sub>` | `main.cpp:115` → `ANALYSIS/analysis.cpp:160-185` | 13 active sub-tools; see `tools-subcommand.md` |

Dispatch is exact `strcmp` (`get_option_cpp.hpp:367`); an unknown command prints a warning and
returns 0 (`:371-373`), so shell scripts do not notice the failure.

## Module dependency sketch

```
main.cpp
 └─ call: NOVA_SV_handler::run (sv_call_main.hpp)
     ├─ RefHandler            reference blocks, 16-mer unique index, RefLocalRepeat (-I), k-mer error profile (-k)
     ├─ BAM_STATUS            SRS insert-size/read-length/depth (-T or computed)
     └─ SV_CALLING_Handler    (SV_core.hpp:1795)
         ├─ somatic branch    -X/-x                       SV_core.hpp:2360-2460
         ├─ FC_BED branch     -b (LRS only)               SV_core.hpp:2464-2492
         └─ basic branch
             ├─ LRS  (-l)     SV_core_LRS.cpp: signals → cluster → MainDBGHandler1 (DBGAssembler.hpp) → Contig_aligner (ksw2) → GT → VCF
             └─ SRS  (-n)     SV_core_SRS*.cpp: SDA/CLP/DRP signals → clusters → dBG/Greedy assembly → ksw2 → GT
                 └─ hybrid    SRS re-genotypes LRS calls; SRS calls near LRS calls are dropped (SV_core_LRS.cpp:2615-2626)
         └─ NovaSVRst.hpp     result items + VCF / PURE_STR writer
```

## Design decisions (as visible in code)

- No threading inside `call`; parallelism is one process per chromosome/region (README.md:88;
  `HitSV_multy_process.sh:192`).
- One reference block at a time, whole chromosome cached as 1 byte/base (`RefHandler.hpp:314-321`).
- Diploid model throughout: at most two haplotype contigs per region (`SV_core_LRS.cpp:1265-1269`).
- Repeat awareness in practice = (a) the SRS-only R file turning R ≥ 40 windows into VNTR candidates
  (`RefHandler.hpp:123`, `SV_core_SRS_ins_del.cpp:663-671`), and (b) depth-scaled thresholds.
  The paper's "empirical signature distributions" survive only as hard-coded constants (see
  `control-flow-*.md`); no distribution file is loaded at runtime.
{% endraw %}
