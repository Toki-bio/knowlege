# CLsat_workflow (satpipe): method profile for comparison with HitSV

Repo: `CLsat_workflow` (Toki-bio/CLsat_workflow), `main` @ `f77448a`, read 2026-10-16, read-only.
Citations are `file:line` from direct reads. Paths under `satpipe/stages/` are abbreviated `S/`.
"unverified" marks anything the code does not show.

---

## 0. One-paragraph summary

Satpipe is an **assembly-only** pipeline. It starts from genome FASTA (or TRF `.dat` files), finds tandem-repeat motifs with TRF, clusters them with vsearch, and uses a **hand-curated** cluster-to-family table to pick seed loci. It then extends each seed outward with windowed `ssearch36` searches against a multi-copy consensus, merges the extended intervals that overlap, re-scans every locus with every family consensus, splits composite loci into segments, and arbitrates the family from bitscores. Orthology across assemblies is called from 600 bp flanks searched **flank-against-flank only**. No stage reads sequencing reads, coverage, alignments or k-mers (§3). The tests check that the code is **equivalent** to the earlier scripts and that the configuration and viewer are well formed. None of them checks **accuracy**: the toy generator writes a planted-truth table that nothing reads.

---

## 1. Method per stage (driver stage IDs from `satpipe/run.py`)

| Driver stage | Script | Input | Tool and command | Key thresholds | Output |
|---|---|---|---|---|---|
| 00_tables | `run.py:219-251` | config | `samtools faidx` | if no `CLsat_cons.rot` is given, the references are derived by min-rotation (`run.py:240-246`) | genomes/families/members.tsv, refs.rot |
| 02_trf | `run.py:277-345` | genome FASTA (or the given `.dat`, `run.py:291-297`) | TRF per scaffold through GNU parallel (`run.py:325-326`) | `trf 2 5 7 80 10 50 2000 -h` (`schema.py:69`; `PIPELINE.md:39`); scaffolds over `large_scaffold_bp` are excluded (`run.py:287-288`) | `<code>.dat`, large_scaffolds.tsv |
| 02b_large ("85mb_recovery") | `S/06_catalogue_finalize/85mb_recovery/run_chunked.sh`, `classify_candidates.py`; `run.py:357-395` | one large scaffold | TRF on 2 Mb chunks with 50 kb overlap (`run_chunked.sh:31-32`); a chunk still running after `TRF_TIMEOUT` is cut into 4 pieces, recursively (`:96-127`); pieces below `MIN_PIECE_BP` are logged as unscanned and never imputed (`:104-107`) | dedup step: same period/copies within 100 bp counts as a duplicate (`:153-155`); candidate kept if period is within ±2 bp of a family mono and copies ≥ 10 (`classify_candidates.py:19-20,60-75`) and its ungapped rotation-aware identity to a consensus is ≥ 0.70 (`:21,37-57,83`) | `<fam>.recovered.bed`, which becomes extra seeds (`run.py:381,389-393`) |
| 03_motifs | `S/01_discovery_clustering/Step3.sh` | `.dat` | awk | period 80–300, copies ≥ 10 (`Step3.sh:19-21,62`) | `<code>.fa`, header `code\|SCF:S-E\|len\|copies\|locuslen` (`:14,63`) |
| 04_filter | `Step4.sh` | motifs | TRF on the motifs themselves (`:29`) | motif dropped as a simple repeat if its internal TRF coverage is > 50 bp (`:21,75`) | `.clean.fas` / `.repeats.fas` |
| 05_rotate | `Step5.sh` | clean motifs | awk, lexicographically minimal rotation | plus strand only by default (`:19,41-42`) | `.rot` |
| 06_cluster | `Step6.sh` | `.rot` | `vsearch --cluster_fast --id 0.90 --strand plus` (`:14-15,37`) | id 0.90 | per-species clusters |
| 07_mafft | `Step7.sh` | clusters with > 1 member (`:27`) | `mafft --localpair --maxiterate 1000` (`:42,108`); rotate to the rightmost start, then realign | singletons skipped, which is why L541 was missed (`S/07_post_catalogue/README.md:20-28`) | `*.rotated.fa` |
| 08_cross | `S/02_seed_expansion_scan/Step8_cross_cluster.sh` | pooled rotated sequences, references prepended | `vsearch --cluster_fast --id 0.90` (`:16,25-27`) | references prepended so they win the centroid election (`:2-6`) | centroids, `cross_clusters/clust_N`; md5 recorded (`run.py:492`) |
| 09_assign | `Step9_family_fasta.sh`; check at `run.py:502-598` | curated `cluster_family.tsv` | `ssearch36` of centroids vs doubled consensuses (`run.py:552`) | **manual** table (`cluster_family.tsv:1-23`); a mismatch is fatal only if the centroid is ≥ 95 % identical at ≥ 100 bits, otherwise a warning (`run.py:578-591`) | `<fam>.members.fa`, `cluster_report.tsv` |
| 10_seeds | `make_bed_from_clsats_headers.sh` | member headers | gawk 3-argument `match` (`:21-25`) | none | seed BED6 (published 750/288/101/44, `PIPELINE.md:103`) |
| 11_expand | `clisat_expand_then_scan_merged.sh` (= `cli_v6.sh`) | seeds, genomes, consensus, MONO | `ssearch36 -n -m 8C -E 1000000` of each window against a 3-copy consensus (`:219-227,294,357`) | step 3×MONO, window 3×MONO, min bits 50, min aligned length 0.67×window, max 500 steps (`:135-138,180-186`); 3 phases: dense walk (`:303-398`), coarse probe of 100×MONO halving down to 50×MONO that **accepts any hit ≥ 50 bits with no length check** (`:409-452`; header `:64-71`), dense pinpoint (`:454-488`); stops on strand flip (`:389-390`), at scaffold ends (`:315-336`) or on signal loss | `expanded_per_locus.bed/tsv` with stop reasons (`:100-109,197`), `dropped.tsv` (`:541-546`) |
| 11 (PASS2) | same script `:548-565` | expanded intervals | `bedtools merge -c 4 -o count,collapse` **with no `-d`** | merge distance 0: only overlapping or touching intervals join | `expanded_merged.bed` |
| 12_master | `S/03_master_arbitration/build_master.sh` | per-family merged BEDs | `bedtools merge -c 4,5 -o distinct,distinct` per species (`:61`), again with no `-d` | merge distance 0 across families | `MASTER.merged.bed` with family provenance |
| 13_rescan | same expansion script `--scan-only` × every family (`run.py:685-`) | MASTER | one `ssearch36` call per window vs the **single** consensus (`cli:593-607`) | window **2×MONO = 292 bp** for 146 bp, step MONO (`cli:138`; `schema.py:105-106`) | `MASTER_<fam>.bits.tsv` |
| 14_arbitrate | `arbitration_patched.sh` | MASTER, bits, seeds | bedtools, python, gawk | DELTA 20, DELTA_AMBIG 15, SEG_BIN 1000, SEG_MIN_RUN 10, SEG_MIN_LOCUS 2000, SEG_SHORT=drop (`:44-47`) | `MASTER.loci.segmented.bed`, `final_assignment.tsv` |
| 15_catalogue | `run.py:723-780` | arbitration + `manual_loci.bed` | scan-only bits for the manual loci | **manual** L541 (`ref/darevskia/manual_loci.bed:1`) | catalogue |
| 16_regions | `build_locus_master.sh`, `S/06_catalogue_finalize/subfamily_rerun/run_extract.sh` | catalogue | samtools | LF/RF 500 bp outside, LB/RB 2×MONO inside (`run_extract.sh:17,63-64`) | `manifest.tsv`, `all_regions.fa` |
| 17_terminal | `terminal_monomers.sh`, `boundary_phase.sh` | regions | `ssearch36` vs doubled consensus | window 1.5 monomers, EDGE_TOL 10 (`terminal_monomers.sh:38-39`) | terminal and phase tables |
| 18_subfamily | `assign_subfamily.sh` | LB window (RB as fallback) | `ssearch36 -E 1000` vs doubled sub-subfamily consensuses (`:41-45,58`) | best bits; confident if Δ ≥ 10 (`build_subfam_term_json.py:138,163`) | `subfamily.json`; empty when `subx/` is missing (`run.py:824-826`) |
| 19_orthology | `S/04_flank_orthology/rebuild_ortho_segmented.py` | catalogue, genomes | `ssearch36` all-vs-all of flanks | see §2.4 | `ortho_pairs.tsv`, `flank_orthologs.json` |
| 20_candidates | `singleton_screen.sh`, `candidate_*.py` (`run.py:853-899`) | Step-6 singletons | `ssearch36 -E 100` vs doubled references (`singleton_screen.sh:61`) | 50-bit floor (`:37`) | C1..Cn, inspection only |
| not in the driver | `edge_refine/refine_edges.py`, `rearbitration/`, `build_consensus_deliverable.py` | — | — | `finalisation.rearbitrate` and `edge_tol` are schema keys (`schema.py:156-157`), but `run.py` never calls `refine_edges` or rearbitration (grep finds no match in `run.py`); `edge_tol` reaches only the viewer (`viewer.py:224`) | — |

**Reads.** No stage takes reads as input. Every input is an assembly FASTA, a TRF `.dat`, or a table derived from them (`schema.py:31-170`).

---

## 2. How loci are defined, refined, assigned and paired across assemblies

### 2.1 Array / locus definition
- A seed is a TRF call with period 80–300 and ≥ 10 copies (`Step3.sh:19-21,62`) whose motif falls in a curated family cluster (`Step9_family_fasta.sh:32-38`).
- Extension runs while each next 3-monomer window scores ≥ 50 bits and aligns ≥ 0.67×window (`cli:384`), stopping at a strand flip (`:389-390`). For arrays longer than about 219 kb, the coarse probe bridges gaps of up to 14.6 kb (`cli:64-71`).
- **Merge distance is 0.** Both merges use `bedtools merge` without `-d` (`cli:560`; `build_master.sh:61`). Gaps are bridged only by the walk itself: one failed window ends Phase 1 (`cli:393-396`).
- **Minimum size.** Only segments of split loci shorter than 2000 bp are dropped (`arbitration_patched.sh:183-186`); an unsplit locus has no length floor. The pipeline has no minimum monomer count beyond the seed's ≥ 10 TRF copies. The viewer alone uses `MIN_MONOMERS_PER_ARRAY = 3` (`satpipe/viewer_html/index.html:1026`; `HANDOFF.md:30`).
- **Locus IDs** are `L<n>` in sorted MASTER order (`arbitration_patched.sh:74`), so they depend on locale and sort (`HANDOFF.md:118`).

### 2.2 Edge refinement
- Driver path: edges sit at the resolution of the walk (3×MONO steps, `cli:180`). PIPELINE says boundaries "differ by ≤200 bp (below the 3×MONO expansion step, i.e. within the method's own resolution)" (`PIPELINE.md:239-241`).
- `refine_edges.py` refines to the base. It builds an end-local consensus from up to 30 terminal monomers with ≥ 0.60 ungapped identity (`:142-176`, `:146`, `:160`), then marks "strict last match" with a run-break at MAXMISS=4 (`:51,246-265`). The result is capped by `ssearch36` window coverage (60 bp windows, 20 bp step; `:211`; `:266`). Scaffold ends within 50 bp are flagged CONTIG_EDGE (`:57,200-201`).
- The thresholds were tuned on hand-checked loci V1, V2 and V517 (`refine_edges.py:8-24`), with no held-out test. The driver does not run this script (§1, last row). Whether the published catalogue used refined edges is **unverified**.

### 2.3 Family and subfamily arbitration (`arbitration_patched.sh`)
- **D3.** Each rescan row becomes a 2×MONO interval (`:101-102`).
- **D3b.** Bins of 1 kb take the family with the highest window score (`:133,149`). A different family that holds ≥ 10 consecutive bins opens a new segment (`:157`). Segments < 2 kb are dropped (published) or merged (`:170-186`).
- **D4/D5.** Per locus, the maximum bits per family (`:238`).
- **D6, three branches:**
  - Branch 1, SEED_CONFLICT or sub-locus: highest bits wins, `_AMBIGUOUS` if winner minus runner-up < 15 (`:275-286`).
  - Branch 2, NO_SEED: same rule (`:288-298`).
  - Branch 3, UNIQUE_SEED: **the seed family is kept even if another family scores higher**, flagged UNIQUE_SEED_CONTRADICTION when Δ ≤ −20 (`:300-311`; `HANDOFF.md:117,30`).
- **Sub-subfamily.** Decided from one 2-monomer window at the array start, RB if LB is missing (`assign_subfamily.sh:36-45`; `run_extract.sh:63-64`), with a confident call at Δ ≥ 10 bits (`build_subfam_term_json.py:138,163`).

### 2.4 Cross-assembly orthology from flanks (`rebuild_ortho_segmented.py`; `clean_ortho_pairs.py` is the reference version)

| Item | Rule | Evidence |
|---|---|---|
| Flank | 600 bp immediately outside each edge, oriented from the array outward | `:53,141-153` |
| Overlap guard | flank discarded (labelled satellite) if it overlaps another catalogued locus; GUARD 0 | `:56,117-122,148-149` |
| Classes | truncated (off the scaffold or < 200 bp); telomeric (≥ 8 copies of TTAGGG/CCCTAA); lowcomplex (distinct 4-mer fraction < 0.35); unique | `:125-133` |
| Residual-satellite drop | flank removed if it hits a doubled **CLsat** consensus at ≥ 60 bp and ≥ 75 % | `:167-183` |
| Search | `ssearch36 -E 1e-3 -b 4000` of the unique flanks **against each other only**, not against the genomes | `:190-191` |
| Pair rule | ≥ 200 bp aligned and ≥ 80 % identity; same-assembly hits excluded; **single best hit per target assembly per query**; no side restriction (L may match R) | `:206-233` |
| Reciprocity | none: an undirected pair is formed from any one directed best hit | `:225-233` |
| Uniqueness in the target genome | not tested | — |
| "Both-sided" | pair supported by an L-sourced AND an R-sourced hit (`viewer_lists.py:119-126`); gives 6 pairs / 11 loci, so the manuscript's "18 of 53" **cannot be reproduced** | `PROVENANCE.md:76-82` |
| Haplotypes | unm and unp have different codes (`configs/darevskia_clsat.yaml:18-19`), so unm↔unp allelic pairs count as "cross-assembly" pairs | `:209` |

There is no PM/MP-style state calling (presence/absence × orthology states); the output is pairs plus a matrix (`:235-251`).

### 2.5 What "85mb_recovery" and "rearbitration" are
- **85mb_recovery.** One 85.83 Mb *unp* scaffold (`JAWDKN010000003.1`) was never scanned because TRF could not finish it (`PIPELINE.md:173-177`). It is now handled by chunked TRF with recursive splitting on timeout (`run_chunked.sh:1-27`) plus motif-identity classification (`classify_candidates.py:1-11`); the driver runs this as stage 02b (`run.py:357-395`).
  - The published run hand-built `sub44_4_4_prefix.fa` (`PIPELINE.md:297`) and left out sub-chunks 44_4_1/2/3 (`PIPELINE.md:255-259`).
  - A ~57 kb telomeric tail was never scanned (`PIPELINE.md:250-253`).
  - The locus reproduces at 422,901 bp against 424,814 bp published (`PIPELINE.md:246-248`).
- **Rearbitration.** Stage 14 re-run on a corrected catalogue (`rearbitration/arbitration_540.sh:1-12`, a wrapper around `arbitration_patched.sh`) after a fresh `--scan-only` pass of every family (`run_pass3.sh:1-33`).
  - The published run reused bits for 439 unchanged loci and recomputed bits for 100 changed ones (`legacy/published_run/06_catalogue_finalize/rearbitration/combine_bits_v2.py:1-5`).
  - The driver does not run it (§1, last row).

---

## 3. Reads, assembly correctness, controls, benchmarks, tests

| Question | Finding | Evidence |
|---|---|---|
| Sequencing reads, coverage, read support, spanning reads | **None.** A grep for minimap/bam/fastq/hifi/depth/kmer/meryl/hifiasm finds no read handling in code. The only "HiFi" strings are in assembly file names | `configs/darevskia_clsat.yaml:18-19` |
| Assembly-correctness checks (collapse, switch, gaps) | None. Assembly N-runs are not handled; only `build_consensus_deliverable.py:122` skips N-containing samples. Scaffold ends are recorded as stop reasons (`cli:102-106`) or CONTIG_EDGE (`refine_edges.py:200-201`), not tested | — |
| Haplotype / phase validation for unm/unp | None. The parental-affinity statistic compares whole-assembly consensus identity unp:nai vs unp:val, from 600 resampled monomers with seed 5 and 300 bootstraps | `build_consensus_deliverable.py:1-14,88-91,173-186` |
| Negative control (shuffled genome, random flanks) | None in code. ssearch36's internal shuffles are only used to estimate bitscore statistics | `HANDOFF.md:50-52` |
| Positive control | Only anecdotal, in a comment: a *Podarcis muralis* control for the IMO-TaqI/Gallotia check, not run by the pipeline | `imotaq_compare.sh:30` (comment), `:221-229` |
| Simulated / planted data | `tests/toy/make_toy.py` plants arrays (12–60 copies, 2–6 % divergence, random phase and strand), one composite CLsat1+CLsat3 array and one orthologous pair, and writes `truth.tsv` (`:4-8,77-117`). **`truth.tsv` is never read**: grep finds no consumer | — |
| Benchmark vs independent truth | None. The "regression" is reproduction of the published catalogue by the same code. Exact items: centroid md5, 750/288/101/44 seeds, 626 segments, span 99.968 %, 8 loci with boundary differences | `PIPELINE.md:197-244` |
| Pytest (`tests/test_schema.py`, `test_check.py`) | Schema defaults equal the published values, path resolution, unknown keys, cross-checks, check-report levels, tool probing. **No tools and no biology** | `test_schema.py:46-140`; `test_check.py:12-106` |
| Toy equivalence | Original scripts at `4d70b56` vs parameterised scripts on the same toy: structural files must be byte-identical; bits columns may differ by **25 % relative tolerance** (the original differs by up to ~15 % between runs) | `run_toy_equivalence.sh:1-11,28-32`; `compare_noise.py:1-15` |
| Driver toy | Driver output equals the harness chain; a second run skips every stage | `run_driver_toy.sh:43-69` |
| Toy test of the cluster table | The table is built from the clusters that contain the reference records, i.e. automatic, not curated | `run_driver_toy.sh:38-40` |
| Molly second case | Runs end to end. Also uses parental-haplotype assemblies, which are **PacBio CLR**; 0 ortholog pairs | `HANDOFF.md:79-86,201-211` |
| Viewer | List rebuild compared with lists embedded in v2.25; headless-Chrome smoke test | `tests/viewer/README.md:3-12` |
| Not covered by the toy | stage 02 TRF on genomes, finalisation, stage-07 scripts (syntax only) | `HANDOFF.md:53-54` |

---

## 4. Reproducibility

| Aspect | State | Evidence |
|---|---|---|
| Config schema | single source of parameters, published defaults, unknown keys are errors | `schema.py:1-7,31-170,272-275`; `README.md:82-84` |
| Hard-coded server paths | the Darevskia config points at `/staging/...` and is marked **unverified**; the therioserver config uses `/home/toki/...` | `configs/darevskia_clsat.yaml:3-9,12-43,103-104`; `configs/darevskia_clsat_therioserver.yaml:1-8` |
| Stale path note | the 07 README still says the scripts carry `/staging` paths, but they are now parameterised | `S/07_post_catalogue/README.md:53-57` |
| `legacy/` | hard-coded paths | e.g. `combine_bits_v2.py:9-11`, `run_pass3_delta.sh:3-6` |
| Inputs not deposited | `CLsat_cons.rot` and the `subx/` sub-subfamily consensuses (therioserver derives refs.rot instead and writes an empty subfamily table) | `PIPELINE.md:296-297`; `configs/darevskia_clsat_therioserver.yaml:4-7` |
| CI | unit job (pytest + `check --no-files --no-tools`) and toy job (conda env, equivalence + driver) | `.github/workflows/ci.yml:1-30` |
| Versions | exact versions only in docs; `environment.yml` pins **minimums** (`trf>=4.09`, `vsearch>=2.22`, …, no version for fasta3/emboss/parallel/gawk); CI uses Python 3.11 against 3.12.11 published | `PIPELINE.md:186-189`; `environment.yml:39-52`; `ci.yml:10` |
| Run manifest | config, tool versions, git commit, DONE fingerprints (sha256 ≤ 64 MB, else size+mtime) | `run.py:42-59,918-934,177-184` |
| Randomness | consensus deliverable seed 5; toy seed 7 | `build_consensus_deliverable.py:37,91,173`; `make_toy.py:23,26` |
| ssearch36 nondeterminism | bitscores vary between runs (up to ~15 % on the toy) and no seed is set. Thresholds at 50 bits, Δ 15/20 and Δ 10 therefore sit inside that noise band; the published data shows 8 split-vs-merged boundary decisions | `HANDOFF.md:50-53`; `PIPELINE.md:239-244` |
| Manual steps | curated cluster table (`PIPELINE.md:261-265`); judgement call for cluster 62 (`cluster_family.tsv:17`); L541 added by hand (`PIPELINE.md:224-233`); hand-made `sub44_4_4_prefix.fa` (`PIPELINE.md:297`); therioserver table re-derived **from the published catalogue** (`configs/..._therioserver.yaml:24-27`; `tests/darevskia/recurate_from_catalogue.py:1-15`) | — |
| Licence | none | `pyproject.toml:64` |

---

## 5. Weaknesses a reviewer could raise about the method

| # | Weakness | Evidence | Suggested remedy |
|---|---|---|---|
| W1 | **Assembly-only evidence.** Array length, copy number and edges come from the assembly sequence alone. Collapsed or expanded arrays and misjoins at flanks cannot be detected | no read input anywhere (§3); `schema.py:31-170` | Map HiFi reads (minimap2 `-x map-hifi`) to unm/unp; report per-locus depth ratio and reads that span from flank to array edge |
| W2 | **Haplotypes treated as parental subgenomes without a test.** unm/unp are separate "genomes"; the affinity test uses whole-assembly consensus identity | `configs/darevskia_clsat.yaml:18-19`; `build_consensus_deliverable.py:12-13,173-186` | Per-contig/per-window parental k-mer affinity (arm/val/nai-specific k-mers); call a hybrid locus only where reads phase it |
| W3 | **Flank "uniqueness" is not genome-wide.** Flanks are filtered only against CLsat consensuses and simple-sequence heuristics, and searched only against other flanks | `rebuild_ortho_segmented.py:125-133,167-183,190-191` | Search each flank against the whole target genome (plus a TE library such as RepeatMasker/Dfam); require a single genome-wide hit with a clear margin over the second |
| W4 | **No reciprocal best hit.** One directed best hit creates a pair; different query flanks can claim the same target | `rebuild_ortho_segmented.py:218-233` | Require reciprocal best hits, ideally from both flanks, plus synteny (neighbouring unique anchors) |
| W5 | **Allelic unm↔unp pairs counted as orthology.** Only same-code hits are excluded | `rebuild_ortho_segmented.py:209`; codes at `configs/darevskia_clsat.yaml:18-19` | Mark pairs within one individual as allelic; analyse them separately with read-based phase support |
| W6 | **"Both-sided" high-confidence set not reproducible** (18 of 53 claimed; the code gives 6 pairs) | `PROVENANCE.md:76-82`; `viewer_lists.py:119-126` | Define the rule in code, re-derive the number, report it with the threshold sweep (70/80/90 already in `:277-289`) |
| W7 | **Thresholds without calibration.** 50 bits, 0.67 aligned fraction, Δ 15/20/10, 80 % / 200 bp, 0.70 identity, SEG_MIN_RUN 10, SEG_MIN_LOCUS 2000 are inherited defaults; refine_edges was tuned on 3 loci | `cli:135`; `arbitration_patched.sh:44-46`; `schema.py:110-126`; `classify_candidates.py:21`; `refine_edges.py:8-24` | ROC on planted truth (`truth.tsv` already exists) and a shuffled-genome null; report sensitivity per threshold |
| W8 | **Thresholds inside the score-noise band.** ssearch36 bits vary by up to ~15 %, and tests allow 25 % | `HANDOFF.md:50-53`; `run_toy_equivalence.sh:28-32` | Fix the ssearch36 statistics (or use raw SW score / identity×length), or repeat scans and report call stability |
| W9 | **Coarse probe bridges up to 14.6 kb** with any ≥ 50-bit hit and no length test, so separate arrays (or an assembly gap) can be joined | `cli:64-71,446-451`; `HANDOFF.md:123` | Require the dense test and read-spanning support across bridged gaps; report bridged gaps explicitly |
| W10 | **Assembly gaps (N runs) not modelled.** A window over Ns simply fails, which ends or splits arrays at gaps | no N handling in `cli`/`refine_edges` (grep); only `build_consensus_deliverable.py:122` | Flag edges within X bp of an N-run or scaffold end as "assembly-limited", like CONTIG_EDGE (`refine_edges.py:200-201`) |
| W11 | **Edge refinement not in the driver**; config keys `edge_tol`/`rearbitrate` exist but are unused | `schema.py:156-157`; `run.py` (no call); `viewer.py:224` | Wire it as a stage or remove the keys; state which edge set the manuscript used |
| W12 | **Seed family kept against a contradicting rescan** (e.g. L186) | `arbitration_patched.sh:300-311`; `HANDOFF.md:117,212` | Report contradictions as unresolved, or arbitrate them on per-window tracks |
| W13 | **Sub-subfamily decided from one 292 bp window** at one edge | `assign_subfamily.sh:36-45`; `run_extract.sh:63-64` | Assign per window along the whole array; report composition and a switch count |
| W14 | **Plus-strand-only clustering drops singletons**; recovery needed a manual screen (L541) | `Step5.sh:6-11`; `Step6.sh:15`; `S/07_post_catalogue/README.md:20-41` | Make `strand: both` the default for new analyses and report the delta; drop the manual locus |
| W15 | **Curated cluster→family table; circular re-curation.** The therioserver run labels its clusters from the published catalogue, then is compared with that same catalogue | `PIPELINE.md:261-265`; `cluster_family.tsv:17`; `recurate_from_catalogue.py:1-15`; `configs/..._therioserver.yaml:24-27` | Assign families by a rule (best-consensus identity/bits) and report agreement with the curated table, rather than deriving the table from the result |
| W16 | **Validation is reproduction, not accuracy.** The regression compares code output with code output; planted truth is ignored | `PIPELINE.md:197-244`; `make_toy.py:87-117` (no consumer) | Add a truth benchmark: planted arrays with known edges, copy number, composites and orthologs; plus an independent-assembly comparison |
| W17 | **Single assembler / one data type per taxon.** The 6 GenBank assemblies are used as given; the molly test also uses haplotype assemblies (CLR) | `configs/darevskia_clsat.yaml:14-20`; `HANDOFF.md:208` | Re-assemble unm/unp HiFi with a second assembler or hifiasm `--hom-cov` settings; compare locus calls across assemblies |
| W18 | **Large-scaffold recovery has hand steps and unscanned sequence** | `PIPELINE.md:246-259,297`; `run_chunked.sh:104-107` | Report unscanned bp per assembly in the manifest (already partly at `run.py:383-388`) and in the paper |
| W19 | **Reference inputs not deposited** (CLsat_cons.rot, subx/) | `PIPELINE.md:296-297` | Deposit them, or make derivation the documented default |
| W20 | **Environment not pinned** to the published versions | `environment.yml:39-52`; `PIPELINE.md:186-189` | Lockfile or container with the exact versions; CI on Python 3.12 |

The genome-size discrepancy (1.4 vs 1.9 Gb) and the narrow taxon focus are not touched by any code here; they are a data and design matter.

---

## 6. Ideas to borrow from HitSV-style approaches, and where they plug in

| Idea (HitSV-style feature) | What it would measure here | Plug-in point |
|---|---|---|
| **Read-spanning support for array ends** | number of HiFi reads anchored in the unique flank (≥ 1–2 kb) that continue into the array past the called edge; edge accepted only with ≥ k spanning reads | new stage after 15_catalogue / 16_regions (uses `manifest.tsv` LF/RF/LB/RB coordinates, `run_extract.sh:17-63`); feeds `refine_edges.py` as an extra bound |
| **Per-locus depth check for collapse/expansion** | array depth ÷ flank depth (or ÷ genome-wide median); ratio ≈ 2 means collapsed haplotypes, < 0.5 means unsupported or false duplication | new stage after 15_catalogue; add columns to `MASTER.final_assignment.tsv` (driver `run.py:723-780`) |
| **Local reassembly of array loci from reads** | gather reads that touch each locus (flank anchors), assemble locally (hifiasm/flye on the subset), call array length/edges/composition on the local contigs, compare with the assembly call | new stage parallel to 11_expand / 13_rescan; the reassembled contig is scanned by the same `cli --scan-only` (`cli:575-619`) so scoring stays comparable; matches the user's "read-first" array-locus prototype |
| **Haplotype-resolved local assembly / phase proof** | partition locus reads by parental k-mers (arm/val/nai), or by heterozygous flank SNPs, then assemble each phase separately; a hybrid or composite locus is claimed only if the reads of one phase carry it | wraps stage 14 output; addresses W2/W5; per-window parental affinity can be added as a viewer track (`satpipe/viewer.py`) |
| **Depth-aware thresholds** | scale min bits / min support by local coverage; mark low-coverage windows "unassessable", not "signal_lost" | `cli` stop logic (`:380-398`) gets a coverage-aware reason; `refine_edges.py` gains a "low_cov" reason next to CONTIG_EDGE (`:200-201`) |
| **Composite TR-MEI style annotation** | annotate TEs (RepeatMasker/Dfam, the user's SINE libraries) inside arrays and in the 600 bp flanks; label CLsat+TE composite arrays and drop TE-derived "unique" flanks | flank classification in `rebuild_ortho_segmented.py:125-183` (add a TE class); D3b segmentation (`arbitration_patched.sh:109-191`) could take TE intervals as extra "families" |
| **Truvari-like benchmark against an independent assembly** | convert locus calls to interval/SV-style records and compare across two assemblies of the same genome (e.g. hifiasm unm vs a second assembler, or arm/val vs a future T2T), reporting precision, recall and edge distance | new `tests/benchmark/`; reuse `bedtools intersect` patterns from `arbitration_patched.sh:209`; numbers go into `manifest.json` |
| **Planted-truth benchmark** | sensitivity, precision, edge error and composite detection against `truth.tsv` | `tests/toy/make_toy.py` already writes `truth.tsv` (`:87-117`); add a scorer to `run_driver_toy.sh` after `:43` and to CI (`ci.yml:24-30`) |
| **Negative-control scramble** | (a) dinucleotide-shuffled or circularly shifted genome through stages 11–19: expected 0 loci and 0 ortholog pairs; (b) random 600 bp non-flank windows through the orthology step: false-pair rate at 80 %/200 bp; (c) shuffled consensus as a decoy family in arbitration | (a) whole driver with a decoy `genomes:` entry; (b) `rebuild_ortho_segmented.py` with a decoy flank set; (c) one extra `families:` entry (`schema.py:48-59`) |
| **Positive controls** | known arrays (e.g. IMO-TaqI / DarSat2 loci, or planted arrays inserted into real flanks) run end to end; recovery rate and edge error | `imotaq_compare.sh` logic as a stage-07 test; planted insertions into real genome segments in `make_toy.py` |
| **Wet-lab validation hook** | export primer-ready flank pairs for the 53 ortholog pairs and for contradicted/composite loci, for PCR/FISH validation | output of 19_orthology (`ortho_pairs.tsv`, `rebuild_ortho_segmented.py:292-296`) |
| **Call stability under score noise** | repeat rescans N times (or use deterministic scores) and report per-locus flip rate of flags and families | wrapper around 13_rescan/14_arbitrate; `compare_noise.py` already measures the noise band |

---

## Key file paths
- `satpipe/run.py` (driver, stages 00–20)
- `satpipe/schema.py` (all parameters)
- `satpipe/stages/02_seed_expansion_scan/clisat_expand_then_scan_merged.sh` (expansion and scan)
- `satpipe/stages/03_master_arbitration/arbitration_patched.sh` (segmentation and arbitration)
- `satpipe/stages/04_flank_orthology/rebuild_ortho_segmented.py` (orthology)
- `satpipe/stages/06_catalogue_finalize/edge_refine/refine_edges.py` (not in the driver)
- `tests/toy/make_toy.py` (`truth.tsv` written, never scored)
- `PIPELINE.md`, `HANDOFF.md`, `PROVENANCE.md`
