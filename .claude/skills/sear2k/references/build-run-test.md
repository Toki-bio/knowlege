# Build, run, test

## Build

Nothing to build. `chmod +x sear` and put `ssearch36` (FASTA36 package;
SINEderella's release check used 36.3.8g from conda env `sinederella`),
`seqkit`, `bedtools`, `samtools` on PATH; `mafft` for `-b`; `minimap2` for
queries >= 1000 bp (`README.md:17-29`).

## Minimal run that exercises the main path

```sh
cd /scratch/run            # splits and intermediates are written to cwd
sear -k Alu.fa genome.fa 0.8 65 50      # -> gen-Alu.bed (BED7), gen-Alu.bnk
sear -k AluY.fa genome.fa 0.8 65 50     # reuses *.2k.part_*.bnk
rm -f *.2k.bnk *.2k.bnk.fai *.2k.part_* *_s.fai   # as step1_search_extract.sh:127 does
```

Exactly the call SINEderella makes per consensus (`step1_search_extract.sh:111`).
Other modes: `sear --best 50 q.fa g.fa` (MAFFT of top 50), `sear -m q.fa g.fa`
(minus rounds), `sear --force-ssearch long.fa g.fa`, `sear --slop 10,20 ...`.

## Test suite

None in this repo. Test material elsewhere that exercises `sear`:

| Test | What it checks | Where |
|---|---|---|
| `tests/test_mask.sh` | synthetic 700 kb genome, 80 planted copies of family A and 60 of B (8 % mutated, both strands); unmasked step 1 must find 80/80 and 60/60, masked 0/80 and 60/60 at planted coordinates (~20 s) | SINEderella `tests/test_mask.sh:17-60`, `docs/MASKING.md` "Usage" |
| `tests/make_toy_satellite.py` | 8 Mb toy genome, 250 bp toy SINE, 200 dispersed copies at 10 % divergence plus planted satellites/arrays/dimers/segmental duplication; truth.tsv | SINEderella `tests/make_toy_satellite.py` |
| `benchmark_sear.sh` | `sear` per query vs `sear_multi`: wall time and BED identity on up to 3 genomes truncated to 100 Mb and 5 queries | SINEderella `benchmark_sear.sh` |
| RepeatMasker comparison | hg38 chr21 and mm39 chr19 recall/precision of the whole pipeline (not committed as a script) | SINEderella `docs/RELEASE_CHECK_2026-10-05.md` |

A ssearch36-free unit test of the awk coordinate mapping would be cheap to add
(feed a hand-written `-m 8C` table) and is the first thing to write before
replacing the aligner.

## Debugging tips

- Run with `-k` and inspect `plus_to_add_1_<part>.bed` (8 columns, window
  coordinates) before the mapping, and `<bank>.all_hits.txt` / `<bank>.all`
  (offsets) after; `cleanup_query` deletes them on exit, so copy them inside a
  `bash -x` session or comment out the trap (`sear:260`).
- `bash -x sear ...` is readable; the script is a flat pipeline.
- ssearch36 stderr is not captured in this repo (SINEderella's step 2 keeps a
  `ssearch.log`); add `2> ssearch.err` on `sear:396` when a part yields nothing.
- Empty BED with a long query: check whether minimap2 mode engaged
  (`Query >= 1000 bp` in the log) and use `--force-ssearch`.
