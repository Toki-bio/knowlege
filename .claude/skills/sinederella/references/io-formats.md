# Run directory and file formats

Line numbers at `90d8ac4`. The pipeline's "data structures" are these files.

## `run_<ts>/` layout

| Path | Written by | Content |
|---|---|---|
| `manifest.txt` | `SINEderella:754-766` | RUN, inputs, THREADS, CHUNK_BP, FLANK, MODE, `VERSION<TAB>tool<TAB>line` per tool, MASK_* when masked |
| `genome.clean.fa` (+`.fai`), `consensuses.clean.fa`, `consensuses.clean.aliases.tsv` | sanitize / canonicalize | `_`->`@U@` headers; alias map of merged RC duplicates |
| `genome.masked.fa`, `mask.bed`, `mask.stats.tsv` | `build_mask` `:693-719` | masked search copy (same headers/lengths) |
| `step1_search_extract.sh`, `step2_asSINEment.sh`, ..., `tools/*.py` | frozen copies `:768-792` | reproducibility |
| `genome.clean_step1/` | step 1 | see below |
| `step2/step2_output*/` | step 2 (+3 writes into it) | see below |
| `results/` | `build_results_links` `:434-464` | symlinks + `report.html`, `alignments/`, `length_variants/`, `consensus_audit/`, `flank_twins.tsv`, `array_flag.tsv` |
| `*.stdout.log`, `*.stderr.log`, `satellites.log`, `*.log` | orchestrator `tee` | |

### `genome.clean_step1/` (`step1_search_extract.sh`)

| File | Format |
|---|---|
| `run_info.txt`, `genome.path`, `consensus.path` | metadata (`:62-76`) |
| `searches/consensuses.clean.part_<NAME>.fa` | one consensus each (`:82`) |
| `searches/gen-consensuses.clean.part_<NAME>.bed` | `sear` BED7: `chrom start end pctid length strand bits` (post-satellite; original `.before_satellites`) |
| `searches/gen-...part_<NAME>.rawhits.tsv` | same 7 cols, every merged hit with pctid >= 65, any length (`sear:341,417,540-544`) |
| `searches/gen-...part_<NAME>.bnk` | copies +-50 bp flank FASTA (unused downstream) |
| `searches/all_hits.labeled.bed` | `chr start end subfam score strand` (`:190-197`) |
| `searches/regions.by_subfam.bed` | `chr start end name[,name...]` (`:207`) |
| `satellites/indication.tsv`, `loci.bed`, `units.fa`, `excluded_hits.bed`, `exclude_regions.bed`, `<q>.kindA.loci.tsv` | satellite stage (`tools/satellite_stage.py:10-20`) |
| `merged_hits.bed` (+ `.before_satellites`) | `chr start end pctid length strand(s)` merged across consensuses (`:218-234`) |
| `extracted.fasta` | `>chr:start-end(strand)` sequences (`:237`) |
| `sampled_30000.fasta`, `subfam_input/{input.fasta,input_NNN.bnk,input_NNN.cons,input.clw,input.msf,input.chunks.tsv,input.clw.al,subfam_log.txt}` | SubFam round (`:247-313`, `SubFam:5-18`) |

### `step2/step2_output/` (`step2_asSINEment.sh`, `step3_postprocess.sh`)

| File | Columns |
|---|---|
| `all_votes.tsv` | `seq best_cons votes summed_bits` (`step2:239-256`) |
| `assignment_full.tsv` | header; `Sequence Subfamily Bitscore Votes Status Threshold RunnerInfo`; Status in {`assigned`, `rejected_low_bitscore`, `no_unanimous`}; Bitscore = **sum over 10 cycles** (`step2:466-478`; audit D17) |
| `assignment_stats.tsv` | `Subfamily Assigned TopN_Bitscore Threshold` (`:459-462`) |
| `assigned.fasta` (`>id|subfam|bits`), `unassigned.fasta`, `subfamilies/<sf>.fasta` | (`:354-437`) |
| `summary.txt`, `ssearch.log`, `seqkit.log` | |
| `self_bits.tsv`, `self_bits_real.tsv` | `subfam bits` (`step3:55-90`) |
| `sim_scores.tsv` | `seq sim_bs self_bs sim_ratio` (`step3:146-157`) |
| `all_sines.bedlike.ALL.tsv` | 12 cols: `chr start0 end subfam_from_extracted assigned_subfam strand best_bitscore leak_flag conflict_flag note sim_bitscore sim_ratio` (`step3:164-179`) |
| `unassigned.tsv` | soft assignments via the search query (`step3:361-414`) |
| `summary.by_subfam.tsv` | per-subfamily counts, LEAK/CONFLICT counts, sim mean/median (`MANUAL.md` 16.8) |
| `plots/<sf>_pctid.tsv`, `*.png` | step 4 |
| `boundary_refinement.tsv` | `subfam population side boundary_bp status last_identity last_elevated_frac bg_identity bg_frac` (`step7:267`) |

### Later outputs

- `results/alignments/<code>_<sf>_{top100,rand100,subfam}.aln.fa`, `manifest.tsv`, `proposals.tsv`, `continuation.tsv` (step 8a / publish; row 1 `_extended`, row 2 original consensus; `[soft]`, `[array]`, `[twin]` marks).
- `results/length_variants/summary.tsv` (`TWO_VERSIONS / UNLINKED_ENDS / SINGLE_MODE / UNRESOLVED`), `results/consensus_audit/summary.tsv` (`MATCH / DIVERGED / SHORTER / UNSTABLE / SKIPPED`), `results/flank_twins.tsv`, `results/array_flag.tsv`.
- `step1c/{cutoffs.tsv,residual.tsv,set_aside.tsv,residual.fasta,subfam_input/}` (`step1c:37-44`).
- `results/report.html` (`step6_report.py`).

## Identifier conventions

- Copy id = `chr:start-end(strand)` from `bedtools getfasta -s`, with `@U@`
  for `_`; incremental step 2 parses it back to BED (`step2:558-564`).
- Consensus names: first header word; composites carry `_P<n>` / `_C<n>`
  suffixes (`tools/vote_eval/retally.py:8`).
- ssearch36 truncates ids at ~60 characters (`tools/length_variants.py:374`),
  cannot open paths > ~120 characters (`tools/satellite_trf_verify.py:222`).
