---
title: HitSV vs SINEderella, sear, SINE_orth_loc, CLsat workflow
type: notes
---

> Type: notes. This comparison mixes verified facts (each with a `file:line` or document
> citation) with judgement about what is worth borrowing and what reviewers are likely to raise.
> Written 2026-10-16.

# HitSV compared with the five tools: what to borrow, and what reviewers can raise

## Pinned versions

| Tool | Repo | Commit read | Digest |
|---|---|---|---|
| HitSV v2.0.4 | hitbc/HitSV | `93880367` | skill `hitsv` |
| SINEderella | Toki-bio/SINEderella | `296fa334` (2026-10-14) | skill `sinederella` (pinned at `90d8ac4`; changes since then are listed below) |
| sear / sear2k | Toki-bio/sear2k | `923ffb6a` | skill `sear2k` |
| SINE_orth_loc v2.0 | Toki-bio/SINE_orth_loc | `53188679` | none (read for this note) |
| CLsat workflow (satpipe) | Toki-bio/CLsat_workflow | `f77448ad` (2026-10-16) | skill `clsat-revision` (project state) |

Citations are paths inside the named repo at these commits. `S/` stands for
`satpipe/stages/` in CLsat_workflow.

## 1. They answer different questions

| | HitSV | SINEderella | sear | SINE_orth_loc | CLsat workflow |
|---|---|---|---|---|---|
| Question | Which SVs does this *individual* carry relative to the reference? | Where are all copies of each SINE family in this assembly, and which subfamily is each? | Where does this consensus occur in this genome? | Is this SINE locus shared or present/absent between *two species*? | Where are the satellite arrays, which family/subfamily, and which are orthologous across assemblies? |
| Primary evidence | read alignments (LRS, SRS, hybrid) | one assembly | one assembly | two assemblies | several assemblies (TRF + ssearch36) |
| Uses reads | yes, only reads | no | no | no | no (`schema.py:31-170`) |
| Unit of output | VCF record (INS/DEL/INV/BND) with genotype | copy with subfamily label | BED hit | aligned locus triplet with SINE/PM/MP state | array locus with family, edges, flank-ortholog pairs |
| Repeat model | generic: local repetitiveness R (SRS only), depth-scaled gates | the bank of consensuses | one consensus | one SINE consensus | 4 family + sub-subfamily consensuses |
| Truth benchmark | GIAB HG002 Q100, 35 HGSVC, 6 T2T species, Truvari | RepeatMasker on chr21/chr19, planted genomes | search-only planted curve (outside the repo) | simulated pair, regression only (`README.md:32-33`) | toy reproduction; `truth.tsv` is written but never scored |
| Controls | PCR validation of CSV loci | none in the pipeline | none | none | none |

**Consequences.**

- HitSV competes with none of the four tools. It adds a layer they all lack: evidence from reads,
  and genotypes within a species.
- The real overlaps are three:
  - SINE_orth_loc's SINE vs empty-site calls and HitSV's MEI insertion and deletion calls.
  - CLsat arrays and HitSV's tandem-repeat (VNTR) expansions, its complex-SV (CSV) loci and its
    arrays of mobile elements embedded in tandem repeats (MEI–TR arrays).
  - SINEderella's annotation step and HitSV's cohort step, which annotates contigs with
    RepeatMasker.
- There is nothing for sear to borrow beyond the benchmarking discipline.

## 2. Where your tools are already stronger than HitSV

Use these points when a reviewer asks "why not just use HitSV?".

- **Sensitivity to diverged elements.**
  - HitSV finds nothing by sequence family. It never searches for an element; it only sees
    alignment breaks. Mobile elements are labelled after the fact, by RepeatMasker on contigs
    (`cohort_csv_analysis.md:237`). Old, fixed and diverged copies, which are the bulk of a SINE
    landscape, are invisible to it by design.
  - ssearch36 recall stays at 0.99–1.00 up to 40 % divergence on planted copies, with 0 false
    positives at E ≤ 1e-3 (`projects/sine-search-strategy.md`, §3–4).
- **Taxon independence.**
  - HitSV's only element-specific logic is a compiled-in human Alu consensus, always on
    (`sv_call_main.hpp:175`, `SV_core_SRS_inv_tra.cpp:822`).
  - Its short-read insertion joining assumes poly-A tails (`SV_core_SRS_combine_contigs.cpp`, ≥ 15 bp
    homopolymer).
  - Your tools take any bank of consensuses.
- **Deep orthology.**
  - HitSV calls against one reference per individual. It cannot say whether a locus is orthologous
    between two species.
  - SINE_orth_loc's three-way alignment (locus 1, locus 2, consensus) does exactly that, and so do
    CLsat's flank-ortholog pairs.
- **Arrays from kilobase units.**
  - HitSV's array test is a regularity test on RepeatMasker hits: DBSCAN `eps=0.5, min_samples=3`,
    length CV ≤ 0.2, spacing SD ≤ 15 bp or CV ≤ 0.2 (`py/mei_tr_array_single_region_detection.py:120,138-143,312`).
    It has no unit-sequence identity check and no chance null.
  - SINEderella's short-piece kind-B check and its chance null (`tools/satellite_screen.py:74-77`)
    are stronger. So is the flank check recorded in skill `tandem-array-unit-check`.
- **HitSV has its own weaknesses** (skill `hitsv`, `pitfalls.md`):
  - every written long-read record is PASS;
  - the v2.0.4 polishing does not reach the output;
  - the cohort tools cap at 24 reference sequences;
  - the documented cohort commands fail as written.
  - Do not present it as a gold standard.

## 3. Per tool: what to borrow, and what reviewers can raise

Severity: **H** = can block acceptance, **M** = needs an added analysis, **L** = a wording or
methods fix.

### 3.1 SINEderella

**Borrow**

| Idea | Why | Where it plugs in | Effort |
|---|---|---|---|
| Benchmark of one species' two assemblies against each other (hg38 vs CHM13, or two haplotypes of one animal), scored Truvari-style with size, position and sequence tolerance. HitSV's multispecies benchmark uses the same Dipcall-on-T2T idea | answers "no independent truth" | new harness next to `tools/coseg_compare/`; output BEDs of step 1–3 | M |
| Full-pipeline planted ROC: 0–40 % divergence × truncation, scored at firm assignment, not at search | the existing curve stops at search; recall is limited by length (`docs/RELEASE_CHECK_2026-10-05.md:30-33`) | extend `projects/sine-search-bench/simulate.py`; new `tests/controls/` | M |
| Negative controls: dinucleotide-shuffled genome, reversed (not complemented) bank, hits per Gb | the pipeline has none | `tests/controls/` | L |
| Local repetitiveness R (largest repeated k in 600 bp, HitSV `fa_stat.hpp:34-101`) reported per locus as a flag | explains low-identity hits inside simple repeats; context for the satellite stage | column added after the `sear` filter (`sear:441-442`) or in `merged_hits.bed` (`step1_search_extract.sh:223-227`) | L (a ~100-line port) |
| Multispecies release panel with frozen banks and pinned versions, numbers per release | HitSV reports six non-human species; reviewers like a fixed panel | release-check document | M |
| **Pairing, not borrowing.** Feed SINEderella subfamily consensuses as the element library to a read-based caller (HitSV or MELT/xTea) and genotype polymorphic insertions across individuals | gives presence/absence within species and recent activity, which SINEderella cannot produce | downstream of `assignment_full.tsv` | M–H |

**Do not borrow:**

- HitSV's Alu logic or its poly-A heuristics. These would reject (TA)n- or (CT)n-tailed SINEs.
- Depth-scaled thresholds. There is no depth in assembly-only mode.

**What reviewers can raise**

| # | Sev | Objection | Evidence | Answer or fix |
|---|---|---|---|---|
| S1 | H | Thresholds are not calibrated: 65 % identity, 0.8 × length, the 0.45 bitscore rule, 10/10 unanimity, LEAK 0.90, the 85 % kind-B cut | `step1_search_extract.sh:116`; `step2_asSINEment.sh:11`; `step3_postprocess.sh:291-303`; `tools/kindB_flank_check.sh:19`; the authors state it themselves at `tools/separation_test.py:18` | planted ROC (above); report the operating point |
| S2 | H | Copy numbers are biased toward full-length copies. AluJb firm recall is 67 %, B1 43 %, and the misses are fragments under 240 bp, not more diverged copies | `docs/RELEASE_CHECK_2026-10-05.md:22,30-33,60-65` | two-tier output (full-length and fragments ≥ 50 bp at an E cut); state the definition in Methods |
| S3 | H | No comparison with current tools: RepeatModeler2, EarlGrey, AnnoSINE_v2 (published in *Mobile DNA*, 2024) | design only (`MANUAL.md:362,1344`) | run AnnoSINE_v2 and RepeatMasker with the same bank on chr21 and one squamate; report overlap |
| S4 | M | A firm label means "closest bank member". With AluY missing, 1,305 AluY copies went to AluSx | `docs/RELEASE_CHECK_2026-10-05.md:34-38` | novelty test: divergence of a copy vs its subfamily's own spread |
| S5 | M | Stochastic `-z 11` scores. Votes reproduce for 99.37 % of copies, but library size moves firm calls from 49.9 % to 72.4 % | `step2_asSINEment.sh:183-191`; `docs/FAMILY_SUBFAMILY_ASSIGNMENT.md:38-44` | deterministic raw SW score with a margin (step 3 already does this, `step3_postprocess.sh:80`); report a margin per copy |
| S6 | M | E-value per 1/20 part with no `-Z`, so significance depends on genome size | `sear:428` (no `-Z`) | set `-Z`, or report bits with a fixed floor |
| S7 | M | Assembly-only: collapsed repeats, copies in gaps and misassemblies are invisible | no read input in the pipeline | read-depth check on a sample of loci (see 3.4); state it as a limitation |
| S8 | M | Subfamily discovery is partly manual, and its automatic parts benchmark poorly (SubFam ARI 0.006–0.27; COSEG better) | `docs/COMPARISON_LOG.md:46-64,93` | publish the curation rules; benchmark the final bank on simulated subfamilies |
| S9 | L | Published Sicista sets predate the strand-aware merge | `docs/STRAND_AWARE_MERGE.md:20` | rerun, or state the version |
| S10 | L | Version pinning: two diverging `sear` copies, fast mode on an unmerged branch, no container | sear2k vs vendored `sear` (67 diff lines); branch `fast-mode` | one canonical `sear`, tagged release, environment file |

### 3.2 sear / sear2k

**Borrow.** Nothing algorithmic; HitSV never searches for a consensus. Two things are worth doing:

- Take benchmarking discipline from HitSV: a recall curve over divergence × length in the repo, not
  only in `projects/sine-search-bench/`, plus the negative controls above.
- Optionally attach the R flag to each hit (as in 3.1).

**What reviewers can raise**

| # | Sev | Objection | Evidence | Answer or fix |
|---|---|---|---|---|
| R1 | M | Upstream sear2k is behind the vendored copy: 2.5 kb vs 3 kb windows, and copies over 500 bp straddling a window edge are cut or lost | sear2k `sear:398,371-375` vs SINEderella `sear:403-404,374-379`; `docs/AUDIT_2026-10-14.md:41` | sync upstream; cite one commit |
| R2 | M | Two scoring regimes in one bank: ssearch36 below 1 kb, minimap2 `map-ont -k 11 -w 3 -c` above, with different identity definitions | `sear:269-318` | report the regime per consensus; give a separate L1 recall (now 52 of 73 full-length on chr21) |
| R3 | M | A locus's identity is the maximum of the merged HSP identities and its length is the longest HSP, not one locus alignment | `sear:428-442` | score with one chained alignment |
| R4 | L | No tests in the sear2k repo; old scripts reverse the meaning of `-m`; the cleanup glob deletes unrelated files | `projects/sear-audit-2026-10.md` findings 1–3, 10 | archive old scripts; private temp dir; a small test |

### 3.3 SINE_orth_loc

How it works, briefly:

1. Copies closer than 300 bp to another copy are set aside, and only the rest get a 300 bp
   **left** flank (`SINE_orth_loc.bash:70-81`).
2. Each flank is mapped to the other genome with `bwa mem` (`:104,111`).
3. Both directions are pooled and clustered by single linkage (`:133-192`).
4. Clusters are split into doubles (2 loci), multis (3–10) and polys (≥ 11) (`:218`).
5. Each cluster is aligned with the consensus by MAFFT, and `ComPair.sh` classifies it:
   - left flank: more than 150 bp, more than 65 identical positions and more than 65 % identity;
   - shared SINE: more than 65 % identity between loci, more than 100 identical bp to each other
     and to the consensus;
   - the rest: SINE, PM or MP from the right-flank test;
   - PM vs MP: whichever locus is more similar to the consensus (`ComPair.sh:78-110`).

**Borrow**

| Idea | Why | Where it plugs in | Effort |
|---|---|---|---|
| **Within-species genotyping to separate fixed insertions from polymorphic ones.** Run HitSV (LRS or hybrid) or another MEI genotyper on a few resequenced individuals per species; genotype each PM/MP/SINE locus | a reviewer's first question for retroposon phylogenetics is incomplete lineage sorting. Allele frequencies let you drop polymorphic markers and support KKSC-style tests | new step after `orth_<sp1>-<sp2>.tsv`; join on locus coordinates | H (data permitting) |
| Read support for empty sites: one read spanning both flanks of an empty site in that species' own reads; for insertions, reads spanning the junctions | turns "absent in the assembly" into "absent in the genome"; the same rule as the CLsat plan ("empty site proven only by one read spanning both flanks") | per-locus check after ComPair (`SINE_orth_loc.bash:382-408`) | M |
| Flank uniqueness score: genome-wide hit count of the 300 bp flank in its own genome, or HitSV-style R ≥ 40 | `bwa mem` placement is accepted with no MAPQ filter; MAPQ is only written into the name (`:105,112`) | after `:105`/`:112`; filter on MAPQ, or on a unique-hit test | L |
| Planted-truth benchmark scored like Truvari: simulate species 2 from species 1 with known insertions, precise deletions, flank decay and paralogous flanks; score PM/MP/SINE precision and recall | the only test is regression equivalence (`README.md:32-33`) | new `tests/` | M |
| Negative control: scrambled target genome; expect about 0 orthologs | already planned in `projects/sq2-distant-orthologous-loci.md`, never run | wrapper run | L |

**What reviewers can raise**

| # | Sev | Objection | Evidence | Answer or fix |
|---|---|---|---|---|
| O1 | H | No accuracy benchmark or negative control; tested only for equivalence between versions | `README.md:32-33` | planted benchmark plus scramble control (above) |
| O2 | H | Placement uses one 300 bp left flank, with no MAPQ, uniqueness or secondary/supplementary filter. Paralogous flanks can create false orthologs, especially when the flank is itself a repeat | `SINE_orth_loc.bash:78-81,104-113` | MAPQ ≥ 20 or a unique-hit test; map the right flank too and require concordance |
| O3 | H | Homoplasy and lineage sorting are not addressed: no target site duplication (TSD) check, no test for precise deletion, no allele frequency. Parallel insertion at one site or an ancestral polymorphism look like shared or PM states | `ComPair.sh` tests flank and SINE identity only (`:78-110`) | report TSDs for shared loci; within-species genotyping (above); state the assumption |
| O4 | M | Systematic exclusions bias the sample: copies within 300 bp of another copy are dropped (`:70-72`), and poly clusters of ≥ 11 loci are written but never processed (`:218`; only doubles and multis are run, `:352,367`) | as cited | report counts excluded per category; at least run polys to a "not assessed" state |
| O5 | M | Fixed cut-offs (65 identical bp, more than 65 % identity, 150 bp, 100 bp) are uncalibrated. The Gekko × Heloderma pilot failed at about 220 Ma | `ComPair.sh:78-96` | calibrate on planted pairs at several divergences; give the divergence limit as a result |
| O6 | M | PM vs MP decided by which locus is closer to the consensus (`SO > ST`), not by which genome has the gap | `ComPair.sh:100-107` | decide from the gap pattern in the alignment; keep similarity only as a tie-break |
| O7 | M | Assembly-only: an "empty site" may be a collapse or misassembly | no read input | read-spanning check (above) |
| O8 | L | One consensus per run, and the input BED comes from upstream (sear/SINEderella), so upstream errors carry through | `SINE_orth_loc.bash:27,54-65` | cite the upstream version and settings |
| O9 | L | Registry merges loci within 60 bp by single linkage (`sine_registry.py:226-232`), which can chain neighbouring insertions | as cited | report a sensitivity to `--tol` |

Reviewers in this area will also name **GPAC** (presence/absence compiler) and the **2-n-way**
suite, which do multi-genome retroposon presence/absence. Compare against at least one of them.

### 3.4 CLsat workflow (satpipe), with Reviewer 1 in mind

This is the tool where HitSV-style ideas matter most. Several of them directly address Reviewer 1's
requests (`clsat-revision`).

**Borrow**

| Idea | Addresses | Where it plugs in | Notes |
|---|---|---|---|
| **Read-spanning support for array ends**: HiFi reads anchored in the unique flank that run past the called edge; accept an edge only with ≥ k such reads | Reviewer 1 (ii); profile W1/W9 | new stage after 16_regions (`S/06_catalogue_finalize/subfamily_rerun/run_extract.sh:17,63-64` gives LF/RF/LB/RB); its result also bounds `refine_edges.py` | the read-first prototype on val218 already does this (39/47 loci have a read-defined end within 5 kb) |
| **Depth ratio per locus**, array vs flank (about 2 means collapsed, under 0.5 means unsupported). HitSV uses depth as a gate (`SV_core_SRS_ins_del.cpp:907-914`); here use it as the measurement | Reviewer 1 (ii): collapses | column in `MASTER.final_assignment.tsv` (`run.py:723-780`) | cheap once reads are mapped |
| **Local haplotype contigs from reads** at array loci: `HitSV call -b loci.bed -B PURE_STR -p HIFI -l reads.bam` | Reviewer 1 (ii, iv): read-proven array structure per haplotype | run on regions under 30 kb; scan the contigs with the same `--scan-only` scorer (`clisat_expand_then_scan_merged.sh:575-619`) | limits: FC_BED needs the region within one 5 Mb block (`SV_core.hpp:889`) and on reference index < 24 (`SV_core.hpp:876-882`), so subset or reorder the reference; the depth cap skips collapsed arrays (`SV_core_LRS.cpp:2320-2321`); arrays longer than the reads are out of reach |
| **Assigning reads to a parent before local assembly**, using arm/val/nai-specific k-mers or flank SNPs | Reviewer 1 (i, iii, iv) | wraps stage 14; affinity track in the viewer (`satpipe/viewer.py`) | the molly method is in skill `amazon-molly-paper` |
| **TE annotation inside arrays and flanks**: HitSV's MEI–TR array idea, applied to CLsat (RepeatMasker class substrings `SINE` / `Satellite`, `FC_joint_calling.hpp:165-180`) | composite arrays; TE-derived "unique" flanks producing false orthologs (profile W3) | flank classes in `S/04_flank_orthology/rebuild_ortho_segmented.py:125-183`; TE intervals as extra "families" in D3b (`S/03_master_arbitration/arbitration_patched.sh:109-191`) | use your SINE library with `-lib` |
| **Comparing array alleles across individuals while ignoring tandem-repeat length changes** (HitSV `tools FC_joint_ana`) | array length polymorphism within species | after per-individual local contigs | 30 kb contig cap (`FC_joint_calling.hpp:702-710`); needs `chr`-style names for its homozygote test |
| **Planted-truth scoring**: `tests/toy/make_toy.py` already writes `truth.tsv` (`:87-117`), but nothing reads it | profile W7/W16 | scorer in `run_driver_toy.sh` after `:43`; CI (`.github/workflows/ci.yml:24-30`) | low effort, high return |
| **Negative controls**: shuffled or circularly shifted genome, random non-flank windows through orthology, shuffled-consensus decoy family | profile W3/W7 | decoy `genomes:` and `families:` entries (`schema.py:48-59`) | — |

**What reviewers can raise.** The full list W1–W20 is in the method profile, summarised here. The
first four map onto Reviewer 1's existing points.

| # | Sev | Objection | Evidence | Answer or fix |
|---|---|---|---|---|
| C1 | H | Array length, copy number and edges come from assemblies alone, with no reads anywhere | no read input (`schema.py:31-170`) | read depth plus spanning reads (above); the assembly-free CLsat amounts already measured from raw reads support this |
| C2 | H | unm/unp treated as parental subgenomes; the affinity test is whole-assembly | `S/06_catalogue_finalize/build_consensus_deliverable.py:1-14,173-186` | per-window read-based affinity. Done for 100 kb blocks (99.0 % / 99.3 %), so report it, and note the swapped labels |
| C3 | H | Flank "uniqueness" is not genome-wide: flanks are searched only against each other (`ssearch36 uflanks vs uflanks`), with one best hit per target, no reciprocity, and unm↔unp allelic pairs counted as orthologs | `S/04_flank_orthology/rebuild_ortho_segmented.py:188-191,206-233`; `configs/darevskia_clsat.yaml:18-19` | search flanks against whole genomes; reciprocal best hits from both sides; label allelic pairs separately |
| C4 | H | The "18 of 53 both-sided" number cannot be reproduced; the code gives 6 pairs over 11 loci | `PROVENANCE.md:76-82`; `satpipe/viewer_lists.py:119-126` | fix the rule in code; re-derive the number before resubmitting |
| C5 | H | Validation is reproduction, not accuracy; the planted truth is never scored | `PIPELINE.md:197-244`; `tests/toy/make_toy.py:87-117` | scorer (above) |
| C6 | M | The coarse probe bridges gaps of up to 14.6 kb with any hit of ≥ 50 bits and no length test, so it can join separate arrays or cross an assembly gap | `clisat_expand_then_scan_merged.sh:62-72,409-452` | require read-spanning support across bridged gaps; list bridged gaps |
| C7 | M | Thresholds (50 bits, Δ 15/20/10) sit inside the ssearch36 score noise band (up to about 15 %; tests allow 25 %) | `HANDOFF.md:50-53`; `tests/toy/run_toy_equivalence.sh:28-32` | deterministic scores, or repeated scans with flip rates |
| C8 | M | The seed family is kept against a contradicting rescan (UNIQUE_SEED_CONTRADICTION) | `S/03_master_arbitration/arbitration_patched.sh:300-311` | report as unresolved |
| C9 | M | Sub-subfamily is called from one 292 bp window at one edge | `S/05_subfamily_terminal_boundary/assign_subfamily.sh:36-45` | per-window composition along the array |
| C10 | M | Manual steps, including circular re-curation: the therioserver cluster table is derived from the published catalogue | `tests/darevskia/recurate_from_catalogue.py:1-15`; `PIPELINE.md:261-265` | rule-based family assignment; report agreement with the curated table |
| C11 | M | No comparison with satellite tools (TRASH, SRF, ULTRA) or assembly-QC tools (NucFreq, Flagger, VerityMap) | — | run TRASH/SRF on one assembly; run NucFreq or Flagger over array loci (already in `RESEARCH_LOCAL_REASSEMBLY.md`) |
| C12 | L | No licence (`pyproject.toml:64`; no LICENSE file); inputs `CLsat_cons.rot` and `subx/` not deposited; minimum-only version pins | `PIPELINE.md:296-297`; `environment.yml:39-52` | add MIT, deposit, lockfile |

## 4. Objections that apply to all four

1. **No independent truth.** Every number is checked against your own runs or against RepeatMasker.
   HitSV's reviewers accepted it because of GIAB, HGSVC and T2T truth sets plus PCR. The cheapest
   equivalent for you:
   - two assemblies of one species (CHM13 vs hg38, or the two haplotype assemblies of one animal);
   - planted truth;
   - a small PCR set for SINE presence/absence and CLsat array junctions.
2. **No negative controls.** A shuffled genome and a reversed bank take an afternoon and pre-empt a
   standard objection.
3. **Thresholds are not calibrated.** Show one ROC per tool and state the operating point.
4. **Assembly-only evidence.** Read support at a subset of loci is enough to quantify how often the
   assembly is wrong. It does not need to cover every locus.
5. **Stochastic ssearch36 statistics (`-z 11`).** Use deterministic raw SW scores where a decision
   depends on a score margin.

## 5. Suggested order

1. CLsat: score `truth.tsv` and fix C4 (the 18/53 number), before anything else in the response
   letter.
2. CLsat: per-locus read depth plus spanning reads (C1); report the 100 kb read-based affinity (C2).
3. CLsat: genome-wide reciprocal flank orthology, with allelic pairs labelled separately (C3).
4. All tools: shuffled-genome and reversed-bank controls.
5. SINE_orth_loc: MAPQ/uniqueness filter, count the exclusions, planted benchmark (O1, O2, O4).
6. SINEderella: planted ROC through firm assignment, plus comparison with AnnoSINE_v2 (S1–S3).
7. Later, data permitting: genotype SINE_orth_loc markers within species to address lineage sorting
   (O3).

## Sources

- HitSV paper: [Nat Biotechnol 2026, doi:10.1038/s41587-026-03332-9](https://www.nature.com/articles/s41587-026-03332-9)
- AnnoSINE_v2: [Mobile DNA 2024](https://link.springer.com/article/10.1186/s13100-024-00331-y)
- GPAC: [PubMed 25261406](https://pubmed.ncbi.nlm.nih.gov/25261406/)
- 2-n-way suite: [Genome Res 2020, 30:1508](https://genome.cshlp.org/content/30/10/1508)
- KKSC statistics for retroposon insertions: [PLoS Comput Biol 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4788455/)
- Method profiles behind this note: skill `hitsv` (references), skills `sinederella`, `sear2k`,
  `clsat-revision`, `tandem-array-unit-check`; `projects/sine-search-strategy.md`,
  `projects/sear-audit-2026-10.md`, `projects/sq2-distant-orthologous-loci.md`.
