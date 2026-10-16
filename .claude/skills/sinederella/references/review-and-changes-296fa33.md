# SINEderella and sear/sear2k: method profile for comparison with HitSV

Checkouts read: SINEderella `main` at `296fa33` (2026-10-14), sear2k `main` at `923ffb6`. Citations are `file:line` in
those checkouts unless they say otherwise. "digest" = skills `sinederella`/`sear2k` (pinned
to `90d8ac4` / `85421e7`); "bench notes" = `projects/sine-search-strategy.md`; "audit notes" =
`projects/sear-audit-2026-10.md`. HitSV features were taken from a summary; see skill `hitsv` for the verified internals (e.g. HitSV does not identify elements by sequence at all).

---

## 1. What changed in SINEderella since `90d8ac4`

35 commits, `90d8ac4..296fa33`, 73 files, +2,908 / −388 (`git diff --stat`). Relevant items:

| Area | Change | Evidence |
|---|---|---|
| **sear windows** | Vendored `sear` now cuts the genome into **3,000 bp windows, step 2,000** (was 2,500). Reason: copies > 500 bp across a window border were lost or cut (3 of 7 lost, 2 cut → 7 of 7 whole). ssearch36 mode now skips queries with len x Slf > 3,000 bp with a warning. | `sear:403-404`, `sear:374-379`; `docs/AUDIT_2026-10-14.md:41` |
| sear2k upstream | **Not synced**: sear2k `sear` still uses 2,500 bp windows and a 2,500 bp guard. The two `sear` copies differ in 67 diff lines. | sear2k `sear:398`, `sear:371-375`; `diff` count run here |
| **minimap2 mode** (queries >= 1 kb) | Rewritten (`38d811a`): query indexed, genome streamed as windows (step max(len, 2 kb), width step + len + 1 kb), `-x map-ont -k 11 -w 3 -c --secondary=no`; real identity and `AS:i`. chr21 L1.3: 4 → 52 of 73 full-length L1s (65 at any length). | `sear:281-318`; audit notes, "Status" table |
| **Strand-aware locus merge** | `tools/merge_hits.py` is the default: overlapping opposite-strand hits become one locus only if they overlap >= 50 % of the shorter one, otherwise two elements with their own strand. Old `bedtools merge` without `-s` (digest pitfall) kept behind `MERGE_STRAND_AWARE=0`. Sicista: 29,295 of 2,177,501 loci (1.3 %) were mixed `+,-`; 600 of 600 resolved pieces now forward vs 284/317 before. Published Sicista results still come from the old merge. | `step1_search_extract.sh:223-227`; `tools/merge_hits.py:13,96`; `docs/STRAND_AWARE_MERGE.md:3-20` |
| **Satellites / arrays** | (a) Kind B check now uses **short pieces** (150 bp flanks, short units <= 400 bp, 150 bp interior piece; ARRAY = median >= 85 in both); old whole-unit check via `SATELLITE_KINDB_METHOD=units`; 2,000-run cap dropped. rsi: 53 s vs 569 s (fast unit check) vs 13 h 24 min (serial); verdicts agree on 2,725 of 2,739 runs, 14 disagreements at the 85 cut. Timema not rerun. (b) Runs checked in **tiles of 8 hits**, edge tiles in halves, only ARRAY tiles excluded (9-unit array + 8 copies: 9 excluded, 0 copies lost; was 17). (c) Kind-B verdicts keyed by (contig, start, end); Dfam names with `/` and `|` matched. | `tools/kindB_flank_check.sh:5-20`; `docs/SATELLITES.md:346-370`; `tools/satellite_stage.py:42,100,173,183-201`; `docs/AUDIT_2026-10-14.md:21,39` |
| **Resume** | `--add` / `--exclude` write `.step1_incomplete` until the rebuild finishes and `--resume` refuses such a run (before, an interrupted `--add` resumed on the source run's copies with exit 0); `extracted.fasta` written to a temp file, `.step1_complete` marker; rebuild re-applies satellite exclusion. | `SINEderella:1075,1133,1255,1307,1380-1383`; `docs/AUDIT_2026-10-14.md:37-38` |
| **Fast mode** | **Not on `main`.** `sear` on main has no `fasta36` / `SEAR_MODE` (grep count 0). It lives only on branch `fast-mode` (`0ae8191`, = `ef1f27f` + merge of `38d811a`): fasta36 round 1 + ssearch36 on the minus bank. `--families-mode fast` on main is something else (search family consensuses, then split to subfamilies, `tools/fast_split.sh`). | `git ls-remote`; `git show FETCH_HEAD:sear` lines 35-39, 69-71; `SINEderella:50-52,188-190,458-460` |
| **Families layer** | `--families`, roll-up of the vote to families (`tools/families.py`, `tools/rollup_families.py`), near-tie rule. hs21 Alu toy: 10,147 subfamily-firm → 10,522 family-firm. | `docs/FAMILIES.md:1-30` |
| **Separation test** | `tools/separation_test.py`: kNN purity + UPGMA label runs vs shuffled-label null and split-in-two control; "No threshold here is calibrated". | `tools/separation_test.py:1-19` |
| **step 2 / 3 correctness** | `sim_ratio` now from raw SW scores (`-m 10 sw_score`) because `-z 11` bits of a tiny library varied 5-10 % and gave ratios > 1; failed ssearch36 cycle now fatal (was `|| true` → silent `no_unanimous`); LEAK after `--add` fixed (0 → 40 tags). | `step3_postprocess.sh:80-83,131-133`; `step2_asSINEment.sh:188-191`; `docs/AUDIT_2026-10-14.md:12-17` |
| **Report** | Fits 400 px, Okabe-Ito colours, manifest `COMMAND / INVOKED_FROM / PIPELINE_COMMIT / SUBFAM_VERSION / ENV` shown; LEAK column was always 0 (fixed). | `SINEderella:361-363`; `docs/AUDIT_2026-10-14.md:19,43` |
| **flankscan fs9 twins** | Body check: a flank-twin pair kept only if bodies are more similar than random copies; Anolis Squam1 82 %, Python Squam3C 96 % of masked-flank twin pairs were hot spots. | `docs/AUDIT_2026-10-14.md:35`; commit `8c78358` |
| **License** | MIT added (`LICENSE:1-3`, `CITATION.cff:12`). Digest's "no license" is obsolete. sear2k: MIT (`README.md:146-148`). | |
| Tests | pytest 97 passed, `run_toy` 13/13 (commit message of `296fa33`; not rerun here). 19 test files, 98 `test_` matches in `tests/*.py` (grep). | |

Still open per the author: flankscan silent empty results when trf/minimap2/dustmasker missing; `--resume` does not run
step 6 / publish; contig names with `|` or `:` in step 2/3 (`docs/AUDIT_2026-10-14.md:51-60`).

---

## 2. Evidence base

### What exists

| Benchmark | Truth | Numbers | Source |
|---|---|---|---|
| hg38 chr21, AluJb/AluSx/AluY (Dfam) | UCSC RepeatMasker, 13,495 Alu; recall = >= 50 bp overlap | firm recall AluJb 66-67 %, AluSx 84-89 %, AluY 90 %; precision 99.8 % with full bank (84.6 % with AluY absent); subfamily agreement 98.9-99.7 %, AluY 89.5 % | `docs/RELEASE_CHECK_2026-10-05.md:14-40` |
| mm39 chr19, B1_Mus1 + B2_Mm2 | RepeatMasker 14,822 B1, 9,412 B2 | firm recall 43.3 % / 33.4 %, precision 99.7 %, subfamily agreement 99.5 / 100 % | `:53-69` |
| Missed-copy analysis | same | missed copies have the same divergence as recalled (AluJb 16.9 vs 16.7 %), median length 140-153 vs 293-299 bp; 97-99 % < 240 bp; 6,607 of 6,608 missed B1 < 240 bp | `:30-33,65` |
| Planted genome, masking | 700 kb, 2 families x 80/60 copies, 8 % subst. | 80/80 + 60/60 found; masked 0/80 + 60/60 | `tests/test_mask.sh` (digest evaluation-data) |
| Planted search-engine test | 8 Mb, 120 bp family, 20 replicates per identity | ssearch36 20/20 at all identities; blastn word 11 4/20 at 85 %, 0/20 <= 75 % | `docs/MASKING.md:17-34` |
| Planted divergence curve (outside the repo) | chr21 background, ~2,800 copies, divergence uniform 0-40 %, 20 % truncated | ssearch36 recall 0.998 (random 300-mer) and 1.000 (dar_squam1), 0.99-1.00 in every 5 % bin to 40 %; **0 fp at E <= 1e-3** | bench notes §3-4 (`projects/sine-search-bench/`) |
| Real genome tool comparison | ssearch36 as reference (D. valentini 1.67 Gb) | fasta36 ktup 4: 99.8 % of ssearch36 loci; one unchunked ssearch36 call keeps 43 % (FASTA36 MAX_BEST cap) | bench notes §5, §7 |
| Hand-curated subfamily sets | owner's group labels (saq 900/9, ccr 800/8, teu 600/6, dmo 500/5) | ARI: COSEG 0.32-0.45; SubFam chunks 0.19-0.27; SubFam+peel 0.09-0.43 | `docs/COMPARISON_LOG.md:46-64` |
| Simulated subfamily trees | `tools/sim/sim_family.py` known labels, 13 scenarios | COSEG ARI 0.28-0.77; SubFam chunks 0.006-0.17; chunk purity 0.45-0.58 | `docs/COMPARISON_LOG.md:66-120` (purity `:93`) |
| Vote reproducibility | two runs, same input | flat outcome identical for 99.37 % of copies; family outcome 100 % | `docs/FAMILY_SUBFAMILY_ASSIGNMENT.md:34-44` |
| Library-size effect on vote | 351-copy part vs padded to 20,000 with shuffled decoys | firm 49.9 % vs 72.4 % | `docs/FAMILY_SUBFAMILY_ASSIGNMENT.md:38-39` |
| Satellites | Vassetzky 2023 counts; genome-wide TRF | windowed TRF recovers 89.7 % of 2,580 TRF loci; kind B short-piece vs unit check 2,725/2,739 agree | `docs/SATELLITES.md:155-213,361-366` (digest) |
| Own earlier runs (rsi, tbr, rle) | earlier SINEderella runs | counts within ±1 % except satellite removals | `docs/RELEASE_CHECK_2026-10-05.md:71-92` |
| Flank twins on real data | random copy-pair body identity | see §1 | `flankscan/tests/real_twins/README.md` |

### What is NOT benchmarked (verified absences or explicit statements)

1. **Assignment accuracy on data with known subfamily labels at genome scale.** The step-2 vote was never scored against
   simulated labels; the simulation (`tools/sim`) feeds only the SubFam/COSEG discovery comparison
   (`docs/COMPARISON_LOG.md:66-120`). Real-data "subfamily agreement" is against RepeatMasker, i.e. another
   consensus-nearest-neighbour method (`docs/RELEASE_CHECK_2026-10-05.md:39-40`), so it is concordance, not accuracy.
2. **False-positive rate on a shuffled / non-SINE genome inside the pipeline.** The bench notes report 0 fp on a 47 Mb
   human background (bench notes §4) but outside the repo, with 5 kb windows, not with `sear`'s 0.8/65 filter; no
   shuffled-genome or reverse-not-complemented-consensus negative control exists in the repo (grep for "shuffl" finds only
   decoy padding, flank shuffles and the separation test null).
3. **Sensitivity vs divergence with the full pipeline filter** (0.8 length + 65 % identity + satellite stage + vote + 0.45).
   The only divergence curve is search-only (bench notes); the release check shows recall is dominated by length
   (`RELEASE_CHECK:30-33`), but the curve for firm assignment vs divergence is not measured.
4. **Threshold calibration.** 65 %, 0.8, 0.45, 0.90 LEAK, 85 % kind B, 10/10 unanimity are fixed constants
   (`step1_search_extract.sh:116`, `step2_asSINEment.sh:11`, `step3_postprocess.sh:291-303`,
   `tools/kindB_flank_check.sh:19`); no ROC or sweep in the repo. The separation test says so explicitly
   (`tools/separation_test.py:18`, `docs/FAMILIES.md:41`).
5. **Run-to-run variance** is quantified once (99.37 %, one genome pair, `FAMILY_SUBFAMILY_ASSIGNMENT.md:43`); counts
   "±1 %" are stated, not distributions (`docs/COMPOSITES.md:88`). No seed control: `-z 11` is not seedable from the
   pipeline (unverified whether FASTA36 exposes a seed flag).
6. **Independent truth.** No comparison against an independent assembly of the same species, a T2T reference, or
   long-read-validated loci; Konkel 2015 Alu loci referenced only (`docs/COSEG_COMPARISON_DESIGN.md:30`, digest).
7. **Competing discovery/annotation tools end-to-end** (RepeatModeler2, EarlGrey, AnnoSINE_v2): design only, no numbers
   (digest evaluation-data, `MANUAL.md:362,1344`).
8. `benchmark_sear.sh` (sear vs sear_multi) has no committed results (digest).
9. sear2k itself: no tests in the repo (repo has only `README.md`, 3 scripts, 1 backup; `ls`).

---

## 3. Weaknesses a reviewer could raise (evidence → one-line remedy)

| # | Weakness | Evidence | Remedy |
|---|---|---|---|
| W1 | **Fixed 65 % identity / 0.8 length rule, uncalibrated, and the length rule is the recall bottleneck**: fragments are invisible, so copy numbers and age spectra are biased to full-length copies. | `step1_search_extract.sh:116` (`sear -k ... 0.8 65 50`); `sear:441-442`; AluJb firm recall 67 %, B1 43 %, misses < 240 bp (`RELEASE_CHECK:22,30-33,60-65`) | Report a two-tier count (full-length >= 0.8 and fragments >= e.g. 50 bp at an E-value cut) and calibrate both cut-offs on a planted ROC. |
| W2 | **Identity is the max of merged HSP identities, length the longest single HSP**, not the locus alignment; a fragmented old copy can fail although it is found. | `sear:428-442` (`-o max,max,distinct,max,max`, `$4>=flen`); bench notes §8.1 (fasta36 lost 18 % for this reason) | Score loci by one chained/global alignment to the consensus (or summed HSP coverage). |
| W3 | **E-value per part, no `-Z`**: significance depends on part size (genome size / 20), so thresholds are not comparable across genomes. | `sear:428` (`-E 2`, no `-Z`; grep `-Z` count 0); flankscan uses `-Z 20000` for this reason (digest aligner-invocation E1) | Fix `-Z` (e.g. genome windows) or report bits with a calibrated floor. |
| W4 | **Window geometry**: 3 kb windows (2.5 kb in sear2k) cap ssearch36 hit length; longer elements need minimap2 mode with a different scorer and identity definition, so one bank mixes two scoring regimes at the 1 kb switch. | `sear:269-276,374-379,403`; sear2k `sear:371-375,398` | Document the regime switch per consensus in output, and sync sear2k upstream to the 3 kb fix. |
| W5 | **Stochastic `-z 11` statistics**: votes are repeated 10x to absorb shuffle noise; outcomes reproduce for 99.37 %, counts move ~1 %; library size changes firmness (49.9 vs 72.4 %). The vote's own noise is used as a filter, not an uncertainty estimate. | `step2_asSINEment.sh:183-191`; `FAMILY_SUBFAMILY_ASSIGNMENT.md:38-44`; `docs/COMPOSITES.md:88` | Use a deterministic score (raw SW, as step 3 already does, `step3_postprocess.sh:80`) with a margin rule, and report a posterior/margin per copy. |
| W6 | **0.45 x 10th-best bitscore rule** is relative per consensus and on summed bits; meaning changes with consensus length, copy number and score scale; it rejects B2 relatives 75-90 % identical. | `step2_asSINEment.sh:11,296-300`; `RELEASE_CHECK:65-69` | Replace with a length-normalised score (sim_ratio) threshold calibrated on simulated subfamilies. |
| W7 | **A firm label means "closest bank member"**: a missing subfamily is absorbed (1,305 AluY → AluSx; ~2,000 AluJ → AluSx firm 10/10). | `RELEASE_CHECK:34-38` | Add an outlier/novelty test (copy-to-consensus divergence vs the subfamily's own distribution) before declaring firm. |
| W8 | **Subfamily discovery is partly manual** and its automatic parts are weak on benchmarks (SubFam chunks ARI 0.006-0.27, purity 0.45-0.58; COSEG better). | `docs/COMPARISON_LOG.md:46-64,93`; digest SKILL fact 8 (`MANUAL.md` 6.1) | State the manual step as such, publish curation rules, and benchmark the final (curated) bank on simulated data. |
| W9 | **Assembly-only**: no reads, so collapsed repeats, missing copies in gaps, chimeric contigs and consensus errors in the assembly are invisible; `assembly_qc.sh` is diagnostic only. | `SINEderella:811-827` (digest pipeline-stages step 6); grep for `.bam`/`fastq` in pipeline code: none | Add an optional read-depth check per locus (copy-number collapse) and an assembly-vs-assembly locus concordance. |
| W10 | **No polymorphism**: a single genome gives presence only; no presence/absence across individuals, no allele frequency, no genotypes. "presence/absence" in docs means module presence inside copies. | `RESEARCH_DIRECTIONS.md:29,65,92-93` | Out of scope for this tool; pair with a read-based MEI genotyper (see §4). |
| W11 | **Comparison mostly to own runs / RepeatMasker**: no independent truth (T2T, second assembly, PCR). | §2 items 1, 6 | Benchmark on CHM13 vs hg38 (same species, two assemblies) and score with Truvari-like matching. |
| W12 | **Satellite kind-B check**: validity rests on an 85 % median cut without a calibration curve; 14 disagreements between old and new methods sit exactly at the cut; uniform-null share test flagged every Alu consensus SAT_B (0 verified). Cost was 85 min of 95 on Timema before the short-piece check; new check not rerun there. `trf` on PATH silently toggles the stage. | `docs/SATELLITES.md:361-368`; `RELEASE_CHECK:46-48`; `step1_search_extract.sh:141-150`; digest pitfalls "Added 2026-10-13" | Calibrate the cut on planted arrays at graded unit divergence; make a missing `trf` an error unless `SKIP_SATELLITES=1`. |
| W13 | **Locus merge semantics**: fixed only on 2026-10-13; published Sicista alignments still use the old non-strand-aware merge. | `docs/STRAND_AWARE_MERGE.md:20` | Rerun affected published sets; note the version in methods. |
| W14 | **Reproducibility / pinning**: no container or lock file; FASTA36 version recorded in `manifest.txt` only; two divergent `sear` copies; fast mode on an unmerged branch; vendored SubFam bumped 1.2.4 → 1.5.0. | `SINEderella:361-363` (manifest); diff of the two `sear`s; `git ls-remote` branches; `docs/AUDIT_2026-10-14.md:24` | Ship an environment file / container, one canonical `sear`, and tag releases cited in papers. |
| W15 | **sear2k as a tool**: repo copy behind; older `sear2k`/`sear2k_bak` have opposite `-m` semantics and hard-coded `-T 24`; cleanup globs delete unrelated files; split reuse is directory-global; long contig names corrupt coordinates (`-w 95` id width). | sear2k `sear2k:190` (`-T 24`); audit notes findings 1-3, 10 | Archive old scripts, use a private temp dir, encode window offsets independent of the id width. |
| W16 | **No negative control** in the pipeline (shuffled genome, reversed consensus, unrelated-taxon genome). | §2 item 2 | Run each bank against a dinucleotide-shuffled genome and a reversed (not complemented) bank; report hits per Gb. |

---

## 4. What to borrow from HitSV, where it plugs in, and what not to borrow

| HitSV idea | Borrow for | Plug-in point | Concrete form |
|---|---|---|---|
| **Repeat-aware thresholds keyed on local repetitiveness** (R = largest repeated k in a 600 bp window) | sear filter and satellite stage: relax/tighten the 65 %/0.8 rule and flag loci in low-complexity or tandem context | after `sear:441-442` (per hit, using the flanking 600 bp) or as a column in `merged_hits.bed` before extraction (`step1_search_extract.sh:223-227`) | Compute R per locus; report it; use it to (a) flag hits embedded in simple repeats (W16), (b) explain low-identity hits. Do not make it a hard filter until calibrated. |
| **DBSCAN on element length and spacing with CV criteria** for MEI–tandem composite arrays | replace/complement satellite kind A (gap <= 100 bp, >= 4 hits) and kind B (`array_order.regular_runs`: >= 5 copies, gap <= 6 kb, ratio 5; wide tier 30 kb/10) | `tools/satellite_screen.py:32-36`, `tools/array_order.py:44-48,51-75` | Cluster consecutive hits by (spacing, hit length); array = cluster with low CV of spacing and length. Gives a period estimate and a calibratable statistic instead of three hand-set tiers; keep SINEderella's chance null (`satellite_screen.py:74-77`) as the significance test, which HitSV's description does not mention. |
| **Truvari-style benchmark against an independent assembly** | evaluate search + vote with an external truth | new harness next to `tools/coseg_compare/`; inputs `results/` BEDs | Run on hg38 and CHM13 (or two haplotype assemblies of one non-model species); lift/align loci; score with size/position tolerance and sequence similarity; report recall/precision per divergence and length bin. Also score RepeatMasker the same way. |
| **Multispecies benchmark** | show behaviour across genome sizes and TE loads | release-check doc | Fixed panel (human, mouse, a bat, a squamate, Timema) with frozen banks and pinned versions; numbers per release. |
| **Negative / positive controls** | W16, W1 | new `tests/controls/` | Negative: shuffled genome, reversed bank. Positive: planted copies over a 0-40 % divergence grid with truncation, scored through the full pipeline (search → firm), extending `projects/sine-search-bench/simulate.py`. |
| **Read-based polymorphic insertion genotyping** (complement, not part) | W10: presence/absence across individuals for SINE copies | downstream of `assignment_full.tsv` (firm copies + subfamily consensuses as the MEI library) | Feed SINEderella's subfamily consensuses as the element library to a read-based caller; genotype reference copies (deletion-type) and non-reference insertions; join by locus. This gives polymorphism and ages without changing SINEderella. |
| **Haplotype-resolved local reassembly** | W9: test whether a locus is collapsed or mis-assembled | optional check on a sample of loci (e.g. flankscan twins, `flankscan/fs9_*`) | For loci with reads available: reassemble ±1 kb, confirm the copy and its flanks; flag collapsed arrays (depth ratio). |
| **Depth-scaled thresholds** | only if reads are added (W9) | same optional check | Depth-normalised support per locus; not applicable to assembly-only mode. |
| **PCR validation** | wet-lab confirmation of selected loci (e.g. polymorphic or array calls) | publication workflow | Small validated set (presence/absence across individuals) as gold standard for W10/W11. |

**Do not borrow:**

- **Hard-coded Alu logic** (element-specific lengths, Alu-specific consensus handling): SINEderella's whole point is
  arbitrary banks across taxa (rsi 21 consensuses incl. composites, Timema 595; digest hit-filtering "Inputs").
- **Poly-A / tail heuristics as a detection criterion**: SINE tails vary by family (A-rich, (TA)n, (CT)n etc.);
  SINEderella already orients simple-repeat tails during bank canonicalisation (digest pipeline-stages step 4) and its
  length-variant test handles end modes (`tools/length_variants.py`, digest). A poly-A rule would reject non-Alu SINEs.
- **Seed-based alignment for detection** (if HitSV relies on read aligners / minimizers for element identification):
  measured recall cliffs at 20-30 % divergence for seed tools (bench notes §4); ssearch36/fasta36 must stay the
  detector (`docs/MASKING.md:17-34`).
- **Depth-scaled thresholds in assembly-only mode**: there is no depth.

---

## Unverified / caveats

- HitSV details are taken from the task description.
- Numbers from docs (`RELEASE_CHECK`, `SATELLITES`, `FAMILY_SUBFAMILY_ASSIGNMENT`, `COMPARISON_LOG`) and the bench notes
  are as recorded by their authors; nothing was rerun here (no ssearch36 in this container per digest).
- Test counts (97 passed, 13/13) are from the `296fa33` commit message.
- Whether FASTA36 `-z 11` can be seeded: unverified.
