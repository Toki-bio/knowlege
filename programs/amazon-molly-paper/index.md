---
title: "amazon-molly-paper — digest"
type: reference
topic: Program digests
generated: true
---

> Type: reference · Claude Code skill `amazon-molly-paper` · [source on GitHub](https://github.com/Toki-bio/knowlege/tree/main/.claude/skills/amazon-molly-paper)

**When to use:** Digest of Ricemeyer, Schaefer et al. 2026 Nature "Gene conversion empowers natural selection in a clonal fish species" (Amazon molly, Poecilia formosa) and its code repo esrice/amazon-molly-paper v0.0.3 - how the two parental haplotypes of a hybrid clonal genome were assembled (trio binning with parental-SPECIES short reads; ancestry-blind Hi-C phasing + contig binning), how phase/crossover breakpoints were verified with long reads, data accessions, and what transfers to the Darevskia unisexualis/armeniaca problem. Use when planning haplotype separation, phase-switch checks, or the molly control in the CLsat revision.

---

{% raw %}
# Amazon molly paper (Nature 2026) — digest for the CLsat revision

Paper: Ricemeyer ES, Schaefer NK, Du K, ... Schartl M, Warren WC. Nature 652:398-404 (9 Apr 2026).
https://doi.org/10.1038/s41586-026-10180-9. Local PDF: `C:\work\CLsat\review\s41586-026-10180-9.pdf`.
Code: https://github.com/esrice/amazon-molly-paper, tag v0.0.3 (Zenodo 10.5281/zenodo.17976428),
local clone `C:\work\CLsat\review\amazon-molly-paper-v0.0.3\` (the browser zip
`amazon-molly-paper-v0.0.3.zip.crdownload` in the same folder is an unfinished download — ignore).
Context for why it matters: skill `clsat-revision`.

## The biology in five lines
- P. formosa = one hybrid cross, P. mexicana mother x P. latipinna father, >=100 ka, near Tampico.
  Diploid, clonal (gynogenesis; achiasmatic meiosis, no crossover, no segregation).
- The two haplotypes are LESS related to each other than each is to its sexual parent species.
- Meselson effect seen (asexual haplotypes diverge faster), but no Muller's ratchet (dN/dS same).
- Gene conversion: 46.1 +/- 7.1 Mb (6.26%) of each fish's genome; 20.4% of genome converted in at
  least one of 19 fish; recurs at same loci; no bias to either haplotype; breakpoints near polyA/T.
  Reverts derived alleles 10.6x more often than spreads them.
- Exactly ONE crossover found genome-wide (homologous breakpoints in both haplotypes, at a TAD edge).
- They checked GC and microsatellite content between haplotypes (no major difference). No satellite
  DNA analysis; repeats only RepeatModeler/RepeatMasker (TETools 1.7). -> satellites are open ground.

## How they got two parental haplotypes (the part we need)
1. **Trio binning with parental SPECIES, not the real parents** (assemblies PforHlat, PforHmex).
   HiFi reads of one female (clonal line WLC 7122) sorted with k=21 k-mers unique to P. latipinna vs
   P. mexicana short reads from an earlier study (Warren 2018), using esrice/trio_binning v1.0.0
   (`find-unique-kmers -k21`, `classify-by-kmers`; kmc 3.2.1 inside). Each bin assembled alone:
   `hifiasm -l0` (v0.16.1-r375 — the same hifiasm version as our D. unisexualis assemblies).
   Script: `analyses/assemblies/formosa_trio_binning/1-assemble-contigs.sh`.
2. **Ancestry-blind check** (PforH1, PforH2): `hifiasm --hg-size 750m --h1/--h2 <Hi-C>` (script says
   v0.19.8-r603). Contigs broken at hand-picked points (`breaks_hap*.tsv`), then binned to a parent
   by `programs/molly/separate_haplotypes.py`: minimap2 -x asm10 hap1 vs hap2, each contig linked to
   its best partner in the other haplotype (graph of connected components), then the whole component
   assigned by total matched bases to latipinna vs mexicana assemblies (mapq >= 30).
   Note: `2-curate.sh` applies `breaks_hap1.tsv` to hap2 too, although `breaks_hap2.tsv` exists —
   looks like a typo in the deposited script; not relevant unless we reuse it.
3. Scaffolding: chromap + YaHS (WarrenLab/hic-scaffolding-nf), Juicebox curation, chromosomes named by
   synteny with Xiphophorus maculatus. Parent species assembled with plain hifiasm 0.16.1 + Hi-C.

## How they verified phase / ancestry switches with reads (reviewer's points ii and iv)
`analyses/recombination/recombination.sh` + `programs/molly/pileup2windows.py` (C version in
`programs/c/`): map HiFi reads of each PARENT species to each formosa haplotype (minimap2 map-hifi),
`samtools mpileup -q20 -Q20 -B`, per position call hom-ref (>80% match), hom-alt, het (top two bases
>90%), skip (<10 bases); per 50 kb window mismatch = (hom-alt + 0.5 het)/called, printed if >=70% of
window called. Where the latipinna and mexicana curves cross = candidate switch. Accepted as real only
if >=2 formosa long reads span it in EACH direction (lat-before/mex-after and mex-before/lat-after),
inspected in IGV. This is exactly the test for "is a switch an assembly error or real".

## Other methods worth knowing
- Short-read phasing: WhatsHap read-backed phasing, then polarity of each phase block set by
  variants fixed between parent species (>=2 supporting, 0 contradicting) — `programs/molly/phase_vcf.py`.
- Gene conversion calls: sites fixed different between the parent species where a formosa fish is
  homozygous; tract = >=3 consecutive such sites, >=100 bp, same direction.
- Ancestral genomes by progressive Cactus 2.4.2 (tree with guppy outgroup GCA_904066995).

## Data (verified 2026-10-15 unless marked)
- BioProjects: P. formosa PRJNA1398335, PRJNA1398336; P. mexicana PRJNA1398337; P. latipinna
  PRJNA1398340; parent Hi-C PRJNA614959 (SRR11430098 mex, SRR11430099 lat).
- Runs: formosa HiFi SRR36680768, formosa Hi-C SRR36680767, mexicana HiFi SRR36680757, latipinna HiFi
  SRR36676996 (from the Step 1 report; SRR36680768 confirmed in SRA).
- 2026 assemblies: I could not find them in NCBI Assembly on 2026-10-15 (they may still be processing).
- Older 2021 haplotype assemblies in NCBI: GCA_013036135.2 (latipinna-derived), GCA_013036085.2
  (mexicana-derived), PacBio CLR. Running TRF on these on therioserver:
  /home/toki/clsat_revision/molly_trf/ (see REVISION_HEADER.md).
- 19 formosa, 12 latipinna, 10 mexicana short-read genomes used for population analyses.

## What transfers to Darevskia
- Trio binning works with reads of the parent SPECIES. For D. unisexualis we have Illumina of
  D. raddei nairensis (17_78, 10_78, 1_78) and D. valentini (5_38, 38_75, 33_75) on KIT
  /data/V/toki/Darevskia_reads/. Caveat: unp sits nearer mixta than valentini in our read tests, so
  valentini may be a poor proxy for the real father; reads with no k-mers of either parent will land in
  "unclassified" — that fraction itself is a measurement.
- The 50 kb window mismatch method + ">=2 spanning reads per direction" rule is a ready, published
  standard for checking unm/unp switch points, and for our array-border loci.
- Gene conversion between haplotypes is large in a 100 ka clone: a satellite shared by both halves
  is not by itself evidence of a phasing error.
{% endraw %}
