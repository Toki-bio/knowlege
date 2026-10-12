---
name: sinederella
description: Internals of SINEderella (Toki-bio), the user's own Bash/Python pipeline for genome-wide SINE discovery: ssearch36 search of consensuses against a fragmented genome (via the vendored `sear`), extraction of copies with flanks, 10-cycle ssearch36 voting for subfamily assignment with a 0.45 bitscore rule, satellite/array screening, length-version and consensus audits, HTML report. Use when reading, modifying, running or benchmarking SINEderella, when its ssearch36 command lines / thresholds / run-directory files need explaining, or when deciding whether to replace or accelerate its ssearch36 steps.
---

# SINEderella internals

> Type: reference (facts verified against the pinned commit) with marked `notes`.

## Pinned source

| Field | Value |
|---|---|
| Repo | https://github.com/Toki-bio/SINEderella |
| Version / tag | `CITATION.cff` version "2026.10", released 2026-10-05 |
| Commit | `90d8ac4` (2026-10-09; full SHA in references/architecture.md) |
| Digested on | 2026-10-11 |
| Language / build | Bash orchestrators + step scripts, Python 3 tools (numpy, matplotlib), one C file (`tools/vendor/kmer_order.c`), JS report helpers; no build system; `pytest.ini` |
| License | not stated in the repo (`CITATION.cff:11` points `license-url` at the repo); treat as unverified |

## What it does in three sentences

Step 1 searches every consensus of a bank, one at a time, against a genome with
`sear` (ssearch36 on 2.5 kb windows, >= 65 % identity, hit >= 0.8 x consensus
length), screens the hits for SINE-derived satellites and tandem arrays, merges
all hits across consensuses into loci and extracts them with `bedtools getfasta`.
Step 2 assigns each extracted copy to the consensus that wins ten independent
`ssearch36 -z 11` runs (stochastic shuffle statistics) and whose summed bitscore
passes 0.45 x the 10th-best copy of that consensus; step 3 adds a single-alignment
similarity ratio, LEAK/CONFLICT flags and soft assignments. Later steps rebuild
consensuses, plot divergence, test length versions, audit consensuses, build
publish alignments and an HTML report; `--add`, `--exclude`, `--resume` and a
masked-search mode reuse an existing run.

## Where to look (read on demand)

| Question | File |
|---|---|
| Repo layout, orchestrator modes, which scripts are live vs historical | `references/architecture.md` |
| Numbered trace of step 1 -> 8b, step 1c, satellites, `--add` | `references/pipeline-stages.md` |
| Every ssearch36 / minimap2 / blastn / mafft call: command, parameters, parsing (Q1) | `references/aligner-invocation.md` |
| Inputs, sizes, chunking; all thresholds from raw hit to firm copy (Q2, Q3) | `references/hit-filtering.md` |
| Where the wall-clock goes, measured timings, parallelism, author notes on speed/sensitivity (Q4) | `references/performance.md` |
| Test sets, planted genomes, ground truth, tool comparisons (ssearch36 vs blastn vs minimap2, COSEG, RepeatMasker) (Q5) | `references/evaluation-data.md` |
| Run directory layout and the key TSV/BED columns | `references/io-formats.md` |
| Gotchas and defaults that surprise | `references/pitfalls.md` |
| Install, run, tests, debugging | `references/build-run-test.md` |

Start with `pipeline-stages.md` when the task is unfamiliar. For the
"replace/accelerate ssearch36" question read `aligner-invocation.md`,
`performance.md`, `evaluation-data.md` in that order. The sibling skill
`sear2k` covers the search tool itself in more depth.

## Ten facts worth knowing before touching the code

1. Two ssearch36 workloads dominate: step 1 (`sear`, one call per consensus x
   20 genome parts, `sear:404`) and step 2 (10 cycles x all consensuses as
   query vs 20,000-copy library parts, `step2_asSINEment.sh:181-189`). On a
   2.1 Gb bat genome with 21 consensuses: step 1 61 min, step 2 12 min
   (`docs/RELEASE_CHECK_2026-10-05.md:75`).
2. The search parameters are `ssearch36 -g -3 -T $T -Q -n -z 11 -E 2 -w 95 -W 70 -m 8C query part 3`
   (`sear:404`) and the vote uses the same flags with `-m 8`
   (`step2_asSINEment.sh:186`); `-z 11` is why votes are repeated 10 times
   (`MANUAL.md:1212`).
3. Search inputs: one consensus per `sear` call (`step1_search_extract.sh:82,111`),
   genome sanitised (`_` -> `@U@`, uppercase; `SINEderella:222-314`) and cut
   into 2500 bp windows / 20 parts by `sear`. Optional N-masking of known
   families (`--mask-run`, `tools/mask_genome.py`, `docs/MASKING.md`).
4. Loci are the union of all consensuses' hits: `bedtools merge -c 4,5,6 -o max,max,distinct`
   **without `-s`** (`step1_search_extract.sh:218-219`), so overlapping +/-
   hits become one locus with strand `+,-`, extracted as `+` (`step8a_extract_alignments.sh:57`).
5. Firm assignment = same consensus wins 10/10 cycles AND summed bitscore >=
   0.45 x the 10th-best summed score of that consensus (`step2_asSINEment.sh:267,307-314,335`).
   Bitscores in `assignment_full.tsv` are sums over 10 cycles
   (`docs/AUDIT_2026-10-05.md` D17).
6. Step 3 adds `sim_ratio` = single ssearch36 bitscore / consensus
   self-score (`step3_postprocess.sh:80,128,146-157`), LEAK when runner-up >=
   0.90 of winner (`:288-300`), CONFLICT from `regions.by_subfam.bed` (`:308-339`).
7. The satellite stage runs inside step 1 before extraction
   (`step1_search_extract.sh:136-146`), reads `sear`'s `*.rawhits.tsv` (hits at
   the identity cut, any length) for kind A and `gen-*.bed` for kind B, and
   removes verified array loci (`tools/satellite_stage.py:4-20`).
8. SubFam (k-mer-ordered chunks of 50 copies from a 30,000-copy sample,
   `step1_search_extract.sh:247-261`) only proposes chunk consensuses; the
   subfamily call is made by a person (`MANUAL.md` 6.1) or the peel in the
   separate SINE-discriminator repo. `docs/COMPARISON_LOG.md` shows COSEG
   beating the SubFam+peel route on hand-curated and simulated sets.
9. The author measured and rejected blastn (loses copies < ~85 % identity),
   blastn-short, minimap2 (`-k 11 -w 5` found nothing) and `ssearch36 -S`
   masking; the time of an ssearch36 search is dominated by re-aligning known
   copies, not by scanning bases (`docs/MASKING.md:19-46,56-60`).
10. Env vars `CHUNK_BP` (30000) and `FLANK` (50) are passed to step 1 as
    `SAMPLE_SIZE` and `BINSIZE` (`SINEderella:112-113,837` vs
    `step1_search_extract.sh:7-8`); the sear flank is hard-coded 50 at
    `step1_search_extract.sh:111`.

## Open questions / unverified

- FASTA36 flag semantics (`-g -3` = gap-extend -3, not "forward strand";
  `-w/-W` display-only): my reading, not run here; `MANUAL.md:1218-1230` reads
  them differently. Results are the same either way (both strands searched).
- Whether `sear`'s untuned `-E 2` (no `-Z`) loses hits on hit-dense parts the
  way flankscan's single big search did (`flankscan/HANDOFF.md:150-153`) is
  untested.
- `benchmark_sear.sh` (sear vs sear_multi) has no committed results.
- Server-side artefacts referenced in docs (`~/tmp/...`, KIT `/data/V/toki/`)
  are not in the repo; numbers quoted from docs are as the author recorded them.

## Change log

- 2026-10-11: initial digest at `90d8ac4`.
