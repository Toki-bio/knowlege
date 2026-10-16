# Using HitSV for SINE polymorphism and satellite arrays in non-human genomes

> Type: notes. This is an interpretation for squamate, talpid and other non-model genomes. It is
> built on the verified facts in the other reference files, but the recommendations are judgement,
> not tested results.

## What carries over unchanged

- **LRS mode has no human training data inside it.** The paper's "empirical signature
  distributions" are reduced to fixed constants. In LRS: >10 events or >35 weighted bases per 200 bp
  (ONT_Q20/HIFI only), clip >1000 bp, support ≥ 0.1 × block mean depth (`control-flow-lrs.md`).
  So species transfer depends on read quality and depth, not on a human model.
- **The R file (`srs_fa_stat`) is pure sequence statistics.** It is the largest repeated k (20–180,
  step 20) inside 600 bp windows (`io-formats.md`). It works for any FASTA. Rebuild it per reference,
  and do it on the alignment FASTA with unchanged record order.
- **The cohort CSV and MEI-TR classifiers key on RepeatMasker class substrings** (`SINE`, `LINE`,
  `Satellite`, `Simple_repeat`...), not on Alu/L1/SVA names (`population-csv-mei-tr.md`). A custom
  RepeatMasker `-lib` with headers like `#SINE/...` and `#Satellite/...` plugs straight in.

## What is human-biased and how to neutralise it

| Issue | Effect on squamate/mole data | Mitigation |
|---|---|---|
| Alu force-calling always on (SRS) | Rare exact matches, but any hit is a fake 313–327 bp INS | Drop SRS `INS` with `IMPRECISE` and SVLEN 313–327, or patch `force_calling_ALU = false` at `sv_call_main.hpp:175` and rebuild |
| Poly-A/T heuristics (`ReadHandler.cpp:23-44`; tail joining ≥15 bp homopolymer in `SV_core_SRS_combine_contigs.cpp`) | Helps SINEs with A-rich tails (many tRNA-derived SINEs). Does nothing for tail-less or (CA)n/microsatellite-tailed SINEs | For SINE insertions from short reads, check the ALT sequences against your consensus library (SINEderella/SubFam) and do not trust the length blindly |
| `-P asm5` default | Over-penalises mismatches when the reference is a different species or population | `-P asm10`; use `asm20` when mapping to a congener's assembly |
| 24-sequence caps in `tools csv_region_500w`, FC_BED and genome-context tools | Silent loss on scaffold assemblies | Reorder FASTA and BAM header so the target scaffolds are indices 0–23. Use `chr`-style names for the `_chr` homozygote test in `FC_joint_ana` |
| Diploid hard-coding | Unisexual/hybrid *Darevskia* are diploid and fine. Triploids are not modelled | For triploid clones, call per haplotype on haplotype-resolved assemblies instead |

## Satellite DNA (e.g. CLsat/DarSat): expect blind spots

- The reference 16-mer index keeps unique k-mers only, so satellite arrays give no long-range
  anchors (`pitfalls.md` B).
- Regions whose depth exceeds 5× (SRS) or `MIN(MAX(5×depth,100),1000)` reads (LRS) are skipped
  silently. Collapsed arrays are exactly these regions.
- Contigs must be spanned by reads, and long-range resolution is limited to 1–50 kb. Megabase-scale
  arrays are out of scope. Use read-based copy-number methods (k-mer counting on raw reads, as in
  the CLsat revision) for array size.
- VNTR-style expansions/contractions inside kb-sized arrays are in scope. Adjacent gaps are not
  merged, so normalise with `bcftools norm` or `truvari collapse` before counting copy changes.
- For a quick look at what HitSV can do at an array locus, run `-b loci.bed -B PURE_STR` on long
  reads. This dumps the haplotype contigs, which can then be dot-plotted against the monomer.

## SINE presence/absence genotyping across individuals

1. Use LRS or hybrid where possible. From short reads alone, large non-MEI insertions are unreliable
   (paper) and the MEI reconstruction relies on poly-A joins.
2. Call per sample with `-m 50 -v 0 -f 0` and the right `-p` / `-P`.
3. Merge across samples externally (SURVIVOR or Jasmine), because `ngs_merging` is unusable.
4. Annotate the INS ALT sequences with your SINE library (RepeatMasker `-lib` or ssearch36 against
   the consensus set). Keep calls where most of the ALT is one SINE plus a TSD.
5. Cross-check against orthology-based presence/absence (sq2 / SINE_orth_loc). Agreement on fixed
   loci validates the caller. HitSV-only polymorphic insertions are the new signal.
6. For SINE-in-satellite composites, adapt `py/mei_tr_array_single_region_detection.py`. Group by
   repeat name (column 10) instead of class keyword, and merge fragmented RepeatMasker hits by ID
   first (`population-csv-mei-tr.md`).

## Open

- Runtime and false-positive rate on a satellite-dense squamate chromosome are unknown. Test on
  one chromosome first.
- `-k` polishing appears to be a no-op (`pitfalls.md` #4). If hybrid accuracy matters, compare
  runs with and without `-k` on one region.
