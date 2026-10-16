---
title: "sinederella: Pipeline stages (control flow)"
type: reference
generated: true
---

> Part of the [sinederella](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/sinederella/references/pipeline-stages.md)

{% raw %}
# Pipeline stages (control flow)

Line numbers at `90d8ac4`. `RUN` = `run_<timestamp>/`.

## Full run: `SINEderella genome.fa bank.fa`

1. Parse args (`SINEderella:130-190`); `THREADS` defaults to `nproc`, `CHUNK_BP=30000`, `FLANK=50` (`:108-113`).
2. Create `RUN`, write `manifest.txt` with tool versions (`:748-766`, `tool_versions` `:319-341`), copy step scripts and `tools/*.py` into `RUN` (`:768-792`).
3. `sanitize_fasta` genome -> `genome.clean.fa`: `_` -> `@U@` in headers, uppercase, 60-col wrap, alphabet check, duplicate names fatal (`:222-314`, `:793-795`). Consensuses -> `consensuses.clean.fa` with first header word only (`:798-800`).
4. `canonicalize_consensus_bank.py`: merge RC duplicates (same-orientation pairs only at >= 98 %, lengths within 90 %), orient simple-repeat tail 3' (`:346-361`, `README.md:52`).
5. `build_mask` if `--mask-*`: `tools/mask_from_run.py` / BEDs -> `mask.bed` -> `tools/mask_genome.py` -> `genome.masked.fa`; `export SEARCH_GENOME` (`:693-719`).
6. `assembly_qc.sh` (diagnostic only, `:811-827`).
7. **Step 1** `step1_search_extract.sh consensuses.clean.fa genome.clean.fa $CHUNK_BP $FLANK` (`:837`; the two numbers land in `SAMPLE_SIZE` and `BINSIZE`, `step1:7-8`):
   1. `samtools faidx` genome if needed (`step1:50-53`); `OUTDIR=genome.clean_step1` (`:57`); `run_info.txt`, `genome.path`, `consensus.path` (`:62-76`).
   2. `seqkit split -i` the bank into `searches/consensuses.clean.part_<NAME>.fa` (`:82-88`, renames seqkit < 0.15 `.id_` files).
   3. For each part file, in `searches/`: `sear -k "$query" "$SEARCH_GENOME" 0.8 65 50` (`:109-123`). `-k` keeps the genome splits between queries; they are deleted once after the loop (`:127`). Output per consensus: `gen-consensuses.clean.part_<NAME>.bed` (BED7), `.bnk`, `.rawhits.tsv`.
   4. Satellite stage unless `SKIP_SATELLITES=1` or no `trf`: `python3 tools/satellite_stage.py --searches searches --genome G --cons C --out satellites --threads T [--no-exclude] --exclude-b ${SATELLITE_EXCLUDE_B:-verified}` (`:135-146`). Rewrites `gen-*.bed` minus excluded hits (originals kept as `.before_satellites`), writes `satellites/{indication.tsv,loci.bed,units.fa,excluded_hits.bed,exclude_regions.bed}` (`tools/satellite_stage.py:10-20`).
   5. `all_hits.labeled.bed` (6 cols: chr start end subfam score strand, subfam from the part file header) and `regions.by_subfam.bed` = `bedtools merge -c 4 -o distinct` (`:161-214`); the comma-joined names are later CONFLICT regions.
   6. `merged_hits.bed` = all BEDs `| bedtools sort | merge -c 4,5,6 -o max,max,distinct` (no `-s`) (`:218-219`); remove intervals overlapping `satellites/exclude_regions.bed` (`:227-234`).
   7. `bedtools getfasta -s -fi $GENOME -bed merged_hits.bed > extracted.fasta` (`:237`; the unmasked genome even in masked mode, `:30-33`); count must equal interval count (`:242-244`).
   8. `seqkit sample -n 30000` -> `sampled_30000.fasta` -> `subfam_input/input.fasta`; `SubFam input.fasta 50` (`:247-261`; wrapper `SubFam:1-54` -> `tools/vendor/SubFam.sh` with `-n 50 -t THREADS -a -K`).
   9. `input.clw` (chunk consensuses) + bank -> `mafft --thread T --threadit 0 --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --preservecase` -> `input.clw.al` for the human subfamily review (`:264-313`, `MANUAL.md` 6.1).
8. **Step 2** `step2_asSINEment.sh consensuses.clean.fa genome.clean_step1/extracted.fasta THREADS step2_output` from `RUN/step2/` (`SINEderella:486-488`):
   1. Clean headers; `seqkit seq -g` the copies (`step2:139-140`); split copies into parts of `SPLIT_SIZE=20000` (`:10,153-174`).
   2. For cycle 1..10, for each part: `ssearch36 -g -3 -T T -Q -n -z 11 -E 2 -w 95 -W 70 -m 8 consensus_clean.fa part >> raw_<cycle>.m8` (`:181-189`); per copy keep the best consensus and the best distinct runner-up (`:193-213`) -> `votes.tsv`.
   3. Aggregate `(copy, consensus)` -> votes, summed bitscore, summed runner-up (`:219-235`); `all_votes.tsv` = best vote per copy (`:239-256`); `unanimous.tsv` = copies with votes == 10 (`:259-282`).
   4. Threshold per consensus = `int(0.45 x bitscore of the min(10, n)-th best unanimous copy + 0.5)` (`:288-317`, `MIN_REL_BITSCORE=0.45` `:11`); status `assigned` / `rejected_low_bitscore` (`:321-345`), runner info tag.
   5. Write `assigned.fasta` (headers `>id|subfam|bits`), `unassigned.fasta`, `subfamilies/<sf>.fasta`, `assignment_stats.tsv`, `assignment_full.tsv` (adds `no_unanimous` rows from `all_votes`), `summary.txt` (`:354-515`).
9. **Step 3** `step3_postprocess.sh RUN THREADS` (`SINEderella:507`): self-bits from `ssearch36 -Q -n -z 11 -E 100 -m 8 cons cons` (`step3:80`); per subfamily `ssearch36 ... cons_single.fa copies.fa` (`:128`), best bits per copy -> `sim_scores.tsv` (`seq sim_bs self_bs ratio`, `:146-157`); 12-column ALL table (`:164-200`); LEAK when `runner_ratio >= 0.90` (`:288-300`); CONFLICT via `bedtools intersect` with `regions.by_subfam.bed` (`:308-339`); soft assignment of unassigned loci from the search query that found them (`:361-414`); `summary.by_subfam.tsv`.
10. `rebuild_consensus_bank.py` (`SINEderella:526`), **Step 4** `step4_plots.sh` (`:532`): per subfamily `ssearch36 -g -3 -Q -n -z 11 -E 100 -T T -m 8 cons copies` for % identity histograms (`step4_plots.sh:139`) and MAFFT profiles.
11. `build_results_links` (`:544`): `results/` symlinks to step outputs.
12. Post-steps (`:546-549`): `tools/array_flag.py`, `flankscan/fs9_run.sh` (flank twins; needs jellyfish > 300 Mb), `tools/length_variants_run.py`, `tools/consensus_audit.py`.
13. `publish_run.sh` when `PUBLISH=1`/`--publish` (needs `SPECIES_CODE`), else `step6_report.sh` -> `results/report.html` (`:558-586`).

## `--add new.fa` (`SINEderella:874-1060`)

1. Find source run (`find_latest_run` `:195-217` or `--run`); sanitise new consensuses; `truly_new` = names not already in the bank (`:897-907`).
2. New `run_add_<ts>/`; symlink genome; merge bank; canonicalize; `cp -r` the source `genome.clean_step1` (`:948-966`).
3. In `searches/`: write `consensuses.clean.part_<NAME>.<ext>` for each new name and run `sear -k qfile GENOME 0.8 65 50` (`:992-1010`); delete splits (`:1013`).
4. `satellite_stage_add` with `--only names` (`:726-742`, `:1018`); `rebuild_from_searches` -> labeled BED, `merged_hits.bed`, `extracted.fasta` (`:594-684`, `:1021`).
5. Incremental step 2 if the source has `assignment_full.tsv`: `--incremental --old-assignment ... --new-beds <union of new BEDs>` (`:1036-1059`): Tier 1 = old firm copies not overlapping new BEDs are carried forward; Tier 2 (rest) gets the 10 cycles (`step2:546-690`). ~75 % carried forward in testing (`MANUAL.md:195`).

## `--exclude NAME` (`:1066-1205`)

Copy step1 output, delete that family's `gen-*.bed` and part file (`:1147-1176`), rebuild, full step 2+3.

## `--resume` (`:1212-1341`)

Refresh scripts in the run dir; skip step 1 if `extracted.fasta` exists; skip
2+3 if `summary.by_subfam.tsv` exists; step 4 if plots exist; always rebuild
`results/` links and run the post-steps.

## Standalone late steps

- **Step 1c** `step1c_deplete.sh RUN [P=25] [BIN=20] [SAMPLE=30000]`: residual = unassigned + assigned copies below the P-th percentile of `sim_ratio` per subfamily, only if >= 0.8 x consensus length (`step1c:98-158`); SubFam bin 20 and MAFFT (`:160-180`). For subfamilies too rare to fill a chunk in the 30k sample (`:12-28`).
- **Step 7** `step7_boundary_refine.sh RUN [FLANK_BASE=50] [STEP_BP=50] [MAX_EXT_BP=1000]`: background = 300 random 70 bp windows, fraction of 2000 pairs above 45 % identity (`:46-125`); per subfamily and side, extend in 50 bp steps until the fraction of elevated pairs <= background + 5 (`confirmed`) or 1000 bp (`undetermined`) (`:204-267`) -> `boundary_refinement.tsv`.
- **Step 8a/8b, publish**: plates (top100 / rand100 / subfam) with `mafft --localpair --maxiterate 1000 --ep 0.123 ... --thread 1`, subfamilies in parallel (`step8a:50-59`); `publish/align_for_publish.sh`, border loop, SINE-discriminator overlay (`PUBLISH_WORKFLOW.md`).

## Threading / parallelism model

Serial pipeline; `THREADS` (env) reaches `sear -T`, the satellite stage
(`--threads`, used as number of single-threaded ssearch36 processes,
`tools/satellite_kindB_verify.py:181-219`), step 2/3/4 `-T`, MAFFT. Step 1's
MAFFT and SubFam use `nproc` regardless (audit D8). Consensuses in step 1 and
cycles/parts in step 2 are sequential loops. Genome-level parallelism only via
`SINEderella_multi --parallel`. Step 8a parallelises over subfamilies with
single-threaded MAFFT (`step8a:53-57`).

## Error handling conventions

`set -euo pipefail` everywhere; orchestrator wraps each step in `set +e` with
`tee` logs (`step1.stdout.log` etc.) and dies on non-zero (`SINEderella:836-846`);
optional stages (satellites, plots, report, post-steps) log a WARNING and
continue; ssearch36 calls end in `|| true` and stderr goes to `ssearch.log`
(`step2:188`).
{% endraw %}
