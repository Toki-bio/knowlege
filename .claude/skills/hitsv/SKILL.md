---
name: hitsv
description: Internals of HitSV (hitbc/HitSV v2.0.4, Nature Biotechnology 2026) — repeat-aware structural-variant caller for long-read, short-read and hybrid BAMs with haplotype-resolved local assembly, plus its cohort complex-SV and MEI-tandem-repeat array tools. Use when reading, running, modifying or debugging HitSV, interpreting its VCF/PURE_STR output, checking which paper claims match the code, or judging how to apply it to SINE/MEI polymorphism and satellite or tandem-repeat arrays in non-human (squamate, mole) genomes.
---

# HitSV internals

> Type: reference (facts verified against the pinned commit) with marked `notes`
> (`references/te-satellite-applicability.md`).

## Pinned source

| Field | Value |
|---|---|
| Repo | https://github.com/hitbc/HitSV (formerly gcSV) |
| Version / tag | v2.0.4 (README changelog 2026-06-24); binary still prints 2.0.0 (`src/clib/desc.h:2`) |
| Commit | `9388036` (full SHA in references/architecture.md), 2026-08-12, single squashed commit |
| Paper | Li, Liu et al., "Structural variant detection across sequencing technologies", Nat Biotechnol 2026, doi:10.1038/s41587-026-03332-9 |
| Digested on | 2026-10-16 |
| Language / build | C++ (+ vendored htslib, ksw2), Eclipse makefile in `Release/`, static `-O3`; Python 3 helpers in `py/` |
| License | GPLv3 |

## What it does in three sentences

`HitSV call` walks the reference in 5 Mb blocks (`RefHandler.hpp:30`) and collects alignment
signals from the BAMs. For long reads these are CIGAR indels ≥ `-m`, dense small-variant windows,
and clips >1000 bp. For short reads they are CIGAR divergence, any soft clip, discordant pairs, and
every window that the precomputed `srs_fa_stat` file marks with local repetitiveness R ≥ 40.
Signals are clustered with depth-scaled thresholds, and each cluster's region ±500 bp is
reassembled with a de Bruijn graph (k=20 for LRS, k=60–120 for SRS) into at most two haplotype
contigs. Those contigs are realigned to the reference with ksw2 (`-P asm5/10/20`), and every
INS/DEL/INV/BND, plus nearby SNPs and indels, is written to a VCF with ratio-rule genotypes. In
hybrid mode the long-read contigs win: short reads only re-genotype them and fill gaps more than
300 bp away. Separate `tools` and `py/` scripts cluster per-sample local contigs across a cohort and
detect MEI–tandem-repeat composite arrays from RepeatMasker and TRF annotation.

## Where to look (read on demand)

| Question | File |
|---|---|
| Code layout, subcommands, who calls what | `references/architecture.md` |
| Structs: RefHandler, RefLocalRepeat, BAM_STATUS, SVE, LRS_SIG, result items | `references/data-structures.md` |
| Inputs, read filters, R stat file algorithm and format, BAM-stat JSON, VCF fields, PURE_STR | `references/io-formats.md` |
| Per-block pipeline branches; LRS signals, thresholds, clustering (paper vs code) | `references/control-flow-lrs.md` |
| SRS signals (SDA/CLP/DRP), R usage, VNTR candidates, Alu force-calling, INV/TRA, hybrid rules | `references/control-flow-srs-hybrid.md` |
| dBG assembly, ksw2 presets, polishing, genotyping and QUAL, somatic mode | `references/assembly-realignment-genotyping.md` |
| Cohort CSV workflow (corrected commands), FC_BED contigs, MEI-TR arrays, contig clustering | `references/population-csv-mei-tr.md` |
| `HitSV tools` sub-tools | `references/tools-subcommand.md` |
| `ngs_merging` / `ngs_fc` (do not use) | `references/merging-forcecalling.md` |
| Gotchas, surprising defaults, bugs | `references/pitfalls.md` |
| Build, demo runs, parallel runs, debugging | `references/build-run-test.md` |
| Applying it to SINEs / satellites in non-human genomes (notes) | `references/te-satellite-applicability.md` |
| What to borrow for SINEderella, sear, SINE_orth_loc, CLsat workflow; reviewer objections (notes) | `projects/hitsv-vs-own-tools.md` (repo root) |

Start with `pitfalls.md` "Read these first" before running it. Start with `architecture.md`
before changing code.

## Ten facts worth knowing before touching the code

1. **The wrapper script is broken for SRS and Hybrid.** `HitSV_multy_process.sh:145` calls
   `srs_bam_stat`; the real command is `ngs_bam_stat` (`src/main.cpp:106`). The bad command exits 0
   and leaves `stat.json` empty (tested).
2. **The `-I` R file is SRS-only and uses chromosome indices.** It is keyed by FASTA record index but
   looked up by BAM tid (`fa_stat.hpp:97`, `SV_core_SRS_ins_del.cpp:664`). It is required with `-n`
   (`sv_call_main.hpp:197`) and ignored in pure LRS mode.
3. **R thresholds.** Windows with R ≥ 40 become forced VNTR insertion candidates
   (`RefHandler.hpp:123`, `SV_core_SRS_ins_del.cpp:663-671`). Repeat-dense genomes therefore
   produce many costly assemblies.
4. **Paper vs code.** Indel ≥ `-m` is 25 bp, not 30 (`sv_call_main.hpp:134`). The LRS clip signal
   needs >1000 bp, not 300 (`SV_core_LRS.cpp:2042`). The SRS read gate is max(3, 0.1×depth) plus a
   score gate (`SV_core.hpp:72-79`). The dBG visit cap is tried at 13, then 9, 5, 1.
5. **Alu force-calling is always on and has no switch.** It uses a human Alu consensus
   (`sv_call_main.hpp:175`, `SV_core_SRS_inv_tra.cpp:822`) and emits a synthetic ~320 bp IMPRECISE
   INS whose length is derived from the position (`SV_core_SRS_ins_del.cpp:1801`).
6. **Hybrid `-k` polishing does not reach the output.** The polished buffer is overwritten by the
   unpolished contig (`SV_core_SRS_ins_del.cpp:3206-3210` vs `:3357,3366`; inferred, not run-tested).
7. **Genotypes come from ratio rules.** LRS uses 6:1 (`SV_core_LRS.cpp:2279-2290`). `GT_cal.hpp`
   is dead code. 0/0 calls are never written, so written LRS records are always PASS; filter on
   QUAL (20–500).
8. **Defaults to override in most uses.** `-v 1` adds SNP/INDEL lines. `-f 1` applies random
   "phasing". `-p ONT_Q20`. `-P asm5` is tuned for human, so use `asm10`/`asm20` for other species.
   `-e` does nothing.
9. **No threads.** Parallelise by `-S/-E` chromosome index, which is a 0-based BAM header index.
   Use `-H` to drop the VCF header and `cat` the parts.
10. **Satellite and collapsed-repeat blind spots.** These come from the unique-only 16-mer reference
    index, the depth caps (SRS >5× mean; LRS > MIN(MAX(5×depth,100),1000) reads), and the 1–50 kb
    long-range limit.

## Open questions / unverified

- Somatic mode with `-X` but no `-x` normal BAM: the code path through an uninitialised handler was
  not run.
- Behaviour with overlapping-pair libraries, where the 99th-percentile insert is ≤ 2× read length:
  the breakpoint model size goes ≤ 0 (`bam_stat.hpp:166`). Not run.
- Whether a very long INS overflows the fixed 1 MB VCF line buffer (`SV_core.hpp:2312-2313`). Not run.
- `-b` with only `-n`: not run.
- No generator for the `-k` 9-mer error profile ships in `src/`. Only `demo/9mer_profile` exists.
- SRS and hybrid modes were not run end to end, because the demo SRS BAM is missing from the repo.
  Their facts are from code reading only.

## Change log

- 2026-10-16: comparison with the user's tools added at `projects/hitsv-vs-own-tools.md`.
- 2026-10-16: initial digest at `9388036`. LRS demo, `srs_fa_stat`, the wrapper-script bug and the
  cohort/MEI-TR tool commands were run on a binary built from source.
