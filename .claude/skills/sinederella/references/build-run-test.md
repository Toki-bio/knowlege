# Build, run, test

## Build / install

No build. Clone, make sure the executable bits are set (`ls -l SINEderella sear SubFam step*.sh`),
and put the directory on PATH (the orchestrator prepends its own dir,
`SINEderella:66`, so `sear`, `SubFam`, step scripts resolve from there).

Dependencies (`README.md:75-86`): `ssearch36` (FASTA36; conda env
`sinederella` on the author's servers has 36.3.8g), `mafft` (7.526), `bedtools`,
`samtools`, `seqkit`, `trf` (satellite stage; skipped without it), `gawk`,
`dustmasker` (flankscan), `jellyfish` (flank twins on genomes > 300 Mb),
Python 3 + numpy + matplotlib. EMBOSS only for the old SubFam.
`tools/vendor/kmer_order.c` is compiled by `tools/vendor/SubFam.sh` as needed
(unverified: I did not read that part).

## Minimal run that exercises the main path

```sh
cd /work/mysp
THREADS=16 SINEderella genome.fa bank.fa            # -> run_<ts>/, ~20 min for 47 Mb x 2 consensuses
THREADS=16 SINEderella --add new_subfam.fa          # -> run_add_<ts>/, incremental step 2
SINEderella --exclude AluJb --run run_<ts>
SINEderella --resume --run run_<ts>
SINEderella --mask-run run_<ts>:VES genome.fa other_bank.fa    # search with VES loci as N
SPECIES_CODE=mysp PUBLISH=1 SINEderella genome.fa bank.fa      # + publish_run.sh
SINEderella_multi full genomes.tsv bank.fa --parallel 2 --threads 8
```

Environment switches: `SKIP_SATELLITES`, `SATELLITE_EXCLUDE[_B]`,
`SKIP_CANONICALIZE`, `CANON_MIN_ID`, `SKIP_REBUILD_CONS`, `SKIP_ASSEMBLY_QC`,
`SKIP_LENGTH_VARIANTS`, `SKIP_CONSENSUS_AUDIT`, `FS9_SKIP`, `SEARCH_GENOME`,
`SINEDERELLA_TOOLS`, `STEP2_TMP_BASE`, `CONSENSUS_AUDIT_JOBS`, `DISC`,
`RAW_ALN_BASE`, `PAGES_INDEX` (`SINEderella:49-57`, `PUBLISH_WORKFLOW.md`).

Running one step by hand:

```sh
cd run_X/genome.clean_step1/searches && sear -k consensuses.clean.part_AluY.fa ../../genome.clean.fa 0.8 65 50
cd run_X/step2 && ../step2_asSINEment.sh ../consensuses.clean.fa ../genome.clean_step1/extracted.fasta 16 step2_output
./step3_postprocess.sh "$PWD" 16
step1c_deplete.sh run_X 25 20 30000 16 0.8
step7_boundary_refine.sh run_X; step8a_extract_alignments.sh run_X mysp
```

## Test suite

| Command | Needs | Time |
|---|---|---|
| `python -m pytest` (from repo root; `pytest.ini`) | mafft, gawk, samtools for the end-to-end ones (skipped otherwise); ssearch36 for `test_satellite_long_paths`, parts of `test_satellite_kindB_verify` | seconds-minutes; 62 passed at the release check |
| `bash tests/test_mask.sh [WORKDIR]` | full toolchain | ~20 s |
| `bash tests/toy/run_toy.sh` | SINEderella + SINE_discriminator checkout (`DISCD`) | seconds |
| `bash flankscan/tests/run_tests.sh 2 3 4 5 6` | trf, ssearch36, gawk, dustmasker | toy, 70 checks |
| `bash tools/vendor/tests/test_kmer_order.sh` | cc, python | |
| `benchmark_sear.sh genomes.tsv bank.fa [100] [5]` | sear, sear_multi | per genome: 2 x 5 queries x 100 Mb |
| `tools/sim/run_sim.sh OUT "1 2 3" 4`, `tools/coseg_compare/run_compare.sh saq` | COSEG at `/data/V/toki/coseg`, SubFam, peel | KIT only |

Best worked examples: `tests/test_mask.sh` (whole step 1 on a planted genome,
reads `all_hits.labeled.bed` for truth comparison), `tests/test_satellite_kindB_verify.py`
(fast vs serial ssearch36 paths on planted arrays), `tests/toy/run_toy.sh`
(plates and publish without a genome).

## Debugging tips

- Logs per step in the run dir: `step1.stderr.log`, `step2/step2.stderr.log`,
  `step2/step2_output/ssearch.log` (ssearch36 stderr of all cycles),
  `satellites.log`, `length_variants.log`, `consensus_audit.log`,
  `flank_twins.log`.
- Hit counts per consensus are printed at the end of step 1
  (`step1:316-327`) and after each `sear` call (`:114-122`).
- To inspect raw search hits keep `gen-*.rawhits.tsv` (identity cut only) next
  to `gen-*.bed`; `gen-*.bed.before_satellites` shows what the satellite stage
  removed; `satellites/excluded_hits.bed` lists it.
- Step 2 temp dir `.step2_assign_XXXX` (under `STEP2_TMP_BASE` or cwd) holds
  `raw_<cycle>.m8`, `votes.tsv`, `vote_summary.tsv`, `unanimous.tsv`,
  `thresholds.tsv`; it is removed on exit (`step2:25-26`), so comment the trap
  to keep it.
- `tools/vote_eval/vote_full.sh CONS COPIES OUT` reproduces the vote keeping
  top-6 scores per copy for offline re-tallying with `retally.py`.
- `SATELLITE_KINDB_SERIAL=1` switches the kind-B check to the reference
  implementation for comparison.
