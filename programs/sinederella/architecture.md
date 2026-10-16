---
title: "sinederella: Architecture"
type: reference
generated: true
---

> Part of the [sinederella](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/sinederella/references/architecture.md)

{% raw %}
# Architecture

Full commit SHA: `90d8ac41748b1dfb5104431948874f8df660a975`
(2026-10-09, "Vendor SubFam 1.2.4 (N never enters the final alignment);
wrapper keeps -a; benchmark_sear.sh: seqkit<0.15 split names"). Shallow clone;
earlier history is only visible through `docs/AUDIT_2026-10-05.md` (139
commits 2026-09-26 to 10-04) and `docs/RELEASE_CHECK_2026-10-05.md`.

Repo URL: https://github.com/Toki-bio/SINEderella. Related repos named in
`README.md:95-99` and docs: SINEdb (site), SINE_consensus (consensus builder),
SubFam (vendored 1.2.4 in `tools/vendor/SubFam.sh`), SINE-discriminator
(`peel_features.py`, `verdict.py`, not vendored), Tal (per-species logs).

## Directory map

| Path | Role | Key files |
|---|---|---|
| `SINEderella` (1368 l) | orchestrator: modes full / `--add` / `--exclude` / `--resume`, masks, sanitising, results links, calls every step | `SINEderella:744-868` mode_full |
| `SINEderella_multi` (659 l) | runs the orchestrator on a TSV of genomes, `--parallel N --threads N` | |
| `sear` (702 l), `sear_multi` (416 l) | the search engine (see skill `sear2k`); multi-query variant (one ssearch36 pass per part, standalone) | `sear:404`, `sear_multi:267` |
| `step1_search_extract.sh` (327 l) | split bank, run `sear` per consensus, satellite stage, labeled/merged BEDs, extraction, 30k sample, SubFam, MAFFT | |
| `step2_asSINEment.sh` (900 l) | 10-cycle vote, thresholds, incremental mode; `asSINEment` (440 l) is the older standalone version | |
| `step3_postprocess.sh` (593 l) | self-bits, sim_ratio, ALL table, LEAK, CONFLICT, soft assignment, summary | |
| `step4_plots.sh` + `plot_subfamily.py` | per-subfamily divergence histograms and profiles (ssearch36 identities + MAFFT) | |
| `step6_report.py` (2803 l) + `report_*.py` | `results/report.html` | |
| `step7_boundary_refine.sh`, `step8a_extract_alignments.sh`, `step8b_publish_report.sh`, `publish/`, `publish_run.sh` | boundary extension test, top100/rand100 plates, publish flow (standalone, not in the auto run) | `README.md:14-21` |
| `step1c_deplete.sh` | residual copies -> SubFam round for rare subfamilies (standalone) | |
| `tools/` | satellite stage (4 files), masks, length variants, arrays, consensus audit, canonicalise bank, `vendor/` (SubFam, sine_consensus, kmer_order.c), evaluation harnesses (`sim/`, `coseg_compare/`, `vote_eval/`, `cpg_div/`) | |
| `flankscan/` | bash/awk flank analysis stages fs0-fs9, own toy tests; only `fs9_run.sh` is wired into the pipeline | `flankscan/HANDOFF.md` |
| `tests/` | pytest unit tests, `test_mask.sh`, `toy/`, `bench/` (MAFFT modes) | `pytest.ini` |
| `docs/`, `MANUAL.md`, `PUBLISH_WORKFLOW.md`, `RESEARCH_DIRECTIONS.md`, ... | design notes, measurements, audits, comparison logs | |

Historical / not called by any entry point (`README.md:66-73`, audit D21):
`asSINEment`, `sear_multi`, `benchmark_sear.sh`, `step4_diagnostic.*`,
`step5*`, `extract_alignments.sh`, `sine_consensus*.sh`,
`sine_pairwise_consensus.sh`, `step1b_*` + `*.js`, `tools/composite_scan.py`,
`import_squamata_run.py`, `workflow.html`.

## Module dependency sketch (full run)

```
SINEderella genome.fa bank.fa
  sanitize_fasta (both)            SINEderella:222-314
  canonicalize_consensus_bank.py   :346-361  (RC-duplicate merge, tail orientation)
  build_mask (optional)            :693-719  -> genome.masked.fa, SEARCH_GENOME
  assembly_qc.sh (diagnostic)      :811-827
  step1_search_extract.sh          :837
     seqkit split -i bank -> searches/*.part_<NAME>.fa          step1:82
     for each: sear -k part.fa $SEARCH_GENOME 0.8 65 50          step1:111   [ssearch36 x 20 parts]
     tools/satellite_stage.py (trf, ssearch36)                   step1:136-146
     all_hits.labeled.bed, regions.by_subfam.bed                 step1:164-213
     merged_hits.bed -> extracted.fasta                          step1:218-237
     seqkit sample 30000 -> SubFam input.fasta 50 -> input.clw   step1:247-261
     mafft L-INS-i of chunk consensuses + bank -> input.clw.al   step1:300-308
  run_step2_step3                  :469-589
     step2_asSINEment.sh bank extracted.fasta THREADS            :488   [ssearch36 x 10 cycles x parts]
     step3_postprocess.sh                                        :507   [ssearch36 self + per-copy]
     rebuild_consensus_bank.py                                   :526
     step4_plots.sh                                              :532   [ssearch36 per subfamily, mafft]
     build_results_links                                         :544
     array_flag, flank_twins (fs9), length_variants, consensus_audit   :546-549
     publish_run.sh (PUBLISH=1) or step6_report.sh               :558-586
```

## Entry points

| Command / mode | Entry | Notes |
|---|---|---|
| `SINEderella GENOME BANK` | `SINEderella:1351-1354` -> `mode_full` `:744` | run dir `run_<timestamp>/` |
| `SINEderella --add new.fa [--run DIR]` | `mode_add` `:874` | searches only new consensuses (`:992-1010`), copies step1 output, incremental step 2 (`:1053-1059`) |
| `SINEderella --exclude NAME` | `mode_exclude` `:1066` | deletes that family's BED, rebuilds, full step 2 |
| `SINEderella --resume [--run DIR]` | `mode_resume` `:1212` | skips completed steps by output presence |
| `--mask-bed BED` / `--mask-run RUN[:names]` | `build_mask` `:693` | full mode only (`:1347-1349`) |
| `SINEderella_multi full|add|exclude genomes.tsv bank.fa --parallel N --threads N` | `SINEderella_multi` | per-genome subprocesses |
| `step1c_deplete.sh RUN [P=25] [BIN=20] [SAMPLE=30000] [T] [MINLEN=0.8]` | standalone | residual -> SubFam |
| `step7_boundary_refine.sh RUN [50] [50] [1000]`, `step8a_extract_alignments.sh RUN CODE`, `publish_run.sh RUN CODE` | standalone / publish | |
| `python -m pytest` | `pytest.ini` | unit tests |

## Design decisions the authors made (and why, if stated)

- **Full Smith-Waterman, never seeds.** `docs/MASKING.md:5,56-60`: blastn word
  11 loses copies below ~85 % identity, blastn-short needs a short-word seed,
  minimap2 has "no SINE-scale sensitivity"; "the class of miss full
  Smith-Waterman was chosen to avoid". Masking known families to N is the
  accepted way to speed up a rescan (`docs/MASKING.md:44-46`).
- **Vote instead of trusting one score.** `-z 11` shuffle statistics move
  bitscores by +-1-3 %; 10/10 unanimity filters borderline copies
  (`MANUAL.md:1212`). Measured: 99.37 % of firm/not-firm outcomes reproduce
  between two runs (`docs/FAMILY_SUBFAMILY_ASSIGNMENT.md:43-44`).
- **Flat bank, one label per copy.** Every consensus is an equal competitor;
  a two-layer family-then-subfamily vote was evaluated and found to change no
  subfamily call but to give firm family labels to the 2-4 % soft copies
  (`docs/FAMILY_SUBFAMILY_ASSIGNMENT.md:1-60`); not implemented.
- **No tandem-repeat or dedup pre-filters before search/SubFam**
  (`MANUAL.md:1344-1360`, `MANUAL.md` 6.1.5): the 0.8-length / 65 % rule makes
  tandem-repeat false positives structurally impossible, and SubFam needs the
  redundant copies.
- **Satellites are screened before assignment** because assignment is the
  costly step and array units are most of such a family's copies
  (`docs/SATELLITES.md:89-99`).
- **Reproducibility:** every iterative MAFFT call uses `--threadit 0`
  (`PUBLISH_WORKFLOW.md:121-123`); scripts are copied into the run dir
  (`SINEderella:768-792`); `manifest.txt` records tool versions (`:319-341`).
- **Modular late steps** (7, 8a, 8b, 1c) run against any finished run and are
  not in the auto sequence (`README.md:23-27`, `RESEARCH_DIRECTIONS.md`).
{% endraw %}
