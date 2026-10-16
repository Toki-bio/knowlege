---
title: "sinederella: Inputs, chunking, and every filter between a raw hit and a firm copy (Q2, Q3)"
type: reference
generated: true
---

> Part of the [sinederella](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/sinederella/references/hit-filtering.md)

{% raw %}
# Inputs, chunking, and every filter between a raw hit and a firm copy (Q2, Q3)

Line numbers at `90d8ac4`.

## Inputs

| Input | What | Size in the author's runs |
|---|---|---|
| Genome FASTA (`GENOME_FASTA`, `.gz` accepted) | sanitised to `genome.clean.fa`: headers `_`->`@U@`, first word kept, uppercase, 60-col, alphabet `ACGTUNRYKMSWBDHV-` (`SINEderella:222-314`) | hg38 chr21 46.7 Mb; mm39 chr19 61.4 Mb; bats 2.1-2.3 Gb; Sicista 3.0 Gb (`docs/RELEASE_CHECK_2026-10-05.md:14,53,73,101`; `docs/SATELLITES.md:325`) |
| Consensus bank FASTA | one record per consensus; canonicalised (RC merge >= 98 %, tail orientation) | 2-4 Alu/B1/B2 of 148-312 bp; rsi bank 21 consensuses incl. 7 composites of 300-600 bp and r9 at 105 bp (`docs/FAMILY_SUBFAMILY_ASSIGNMENT.md:99-100`); Mammalia bank 58 (`docs/SATELLITES.md:325`); Timema 595 (`SubFam:29`) |
| Optional mask BED / earlier run | N-masked copy of the genome for the search only (`docs/MASKING.md`) | tbr VES mask 640,519 loci, 136 Mb = 6 % (`docs/MASKING.md:122`) |

### Chunking

- **Genome:** `sear` cuts the search genome into 2500 bp windows with step 2000 (`seqkit sliding -g -s 2000 -W 2500`) and 20 equal-record part files (`seqkit split2 -p 20`) (`sear:379-380`); the split is made once per directory and reused by every consensus (`-k`, `step1_search_extract.sh:107-111`), deleted at `step1:127`. Each ssearch36 call sees one part (~1/20 of windows: for a 2 Gb genome ~50,000 windows, 125 Mb).
- **Query:** one consensus per `sear` call (`step1:82,109-123`); `sear_multi` (standalone) sends the whole bank per part.
- **Copies for the vote:** `extracted.fasta` split into parts of 20,000 records (`step2:10,153-174`); consensuses are the query, copies the library. Library-size effects: the small last part voted 96.0 % unanimous vs 97.2-98.4 % in full parts; padding with shuffled decoys changed 1.3-1.7 % of outcomes (`docs/FAMILY_SUBFAMILY_ASSIGNMENT.md:60-63`).
- **Satellite unit check:** `--max-verify-runs 2000` longest runs per consensus (`tools/satellite_stage.py:89`).

## Filter chain, in order

| Stage | Rule | Code | Effect / evidence |
|---|---|---|---|
| S1 ssearch36 | `-E 2` per part, both strands, `-z 11` | `sear:404` | all local alignments above E 2 |
| S2 window merge | `bedtools merge -s` of hits on one strand within one window; per merged interval keep max hit length, max query span, max %id, max bits | `sear:413` | several local alignments of one copy become one interval |
| S3 identity | merged `%id >= 65` (Homology, positional 4; `step1:111` passes 65) | `sear:417-418` | `rawhits.tsv` keeps everything passing this, any length (satellite stage input) |
| S4 length | merged interval's max hit length `>= int(0.8 x consensus length)` (Slf 0.8, `step1:111`) | `sear:418` | **the dominant recall limit**: 97-99 % of missed RepeatMasker Alu on chr21 are < 240 bp fragments with the same divergence as recalled copies (`docs/RELEASE_CHECK_2026-10-05.md:30-32`); 6,607/6,608 missed B1 are < 240 bp (`:62-66`). Lower to 0.2 for satellite monomers (`docs/SATELLITES.md:163`). |
| S5 genome merge | window -> genome coordinates, `bedtools slop -l 1 -r -1`, `merge -s -c 4,5,6,7 -o max,max,distinct,max` | `sear:521-538` | boundary-straddling copies from two windows collapse; per-consensus BED7 `gen-<q>.bed` |
| S6 satellite exclusion | kind A: TRF loci in +-1 kb windows around raw hits whose unit aligns >= 45 bp to the consensus at E <= 0.01 (`satellite_trf_verify.py:225`); kind B: regularly spaced runs (>= 5 hits, gaps <= 6 kb, within 5x of median gap) with a chance null, `verified` when lag-1/lag-2 unit identity median >= 85 % (`satellite_kindB_verify.py:16`); hits in such loci removed from `gen-*.bed`, regions to `exclude_regions.bed` | `tools/satellite_stage.py:4-20`, `step1:136-146,227-234` | rsi MEG-RS 1,715 -> 129 copies (21 verified arrays, 2,243 hits excluded) (`docs/RELEASE_CHECK_2026-10-05.md:87-90`) |
| S7 locus merge | all consensuses' BEDs: `bedtools sort | merge -c 4,5,6 -o max,max,distinct` (**not strand-aware**) | `step1:218-219` | one locus per genomic interval whatever consensus found it; strand may be `+,-`; `regions.by_subfam.bed` records which consensuses overlap (CONFLICT later) |
| S8 extraction | `bedtools getfasta -s` on the unmasked genome, no flank (the `50` passed to `sear` only affects `.bnk`, which is not used downstream) | `step1:237` | `extracted.fasta` ids `chr:start-end(strand)` |
| S9 vote | 10 cycles; per cycle best consensus by bits; firm candidate iff the same consensus wins 10/10 | `step2:181-282` | ~2-4 % of copies are not unanimous (`docs/FAMILY_SUBFAMILY_ASSIGNMENT.md:29-35`); 99.37 % reproducible |
| S10 bitscore threshold | summed-over-10-cycles bits `>= int(0.45 x bits of the min(10,n)-th best unanimous copy of that consensus + 0.5)` else `rejected_low_bitscore` | `step2:11,288-345` | designed against truncated/degraded copies (`MANUAL.md:1240-1246`); rejects ~2 % on rsi (1,782 of 66,567, `docs/FAMILY_SUBFAMILY_ASSIGNMENT.md:56-58`); B2 subfamilies 75-90 % identical to the bank member fail it (`RELEASE_CHECK:66-69`) |
| S11 soft calls | `no_unanimous` and `rejected_low_bitscore` copies keep their best query as a soft label (`unassigned.tsv`, via the original search query BED) but are not counted as members | `step3:361-414`, `docs/FAMILY_SUBFAMILY_ASSIGNMENT.md:17-19` | |
| S12 LEAK | runner-up summed bits / winner >= 0.90 | `step3:288-300` | flags near-ties between consensuses |
| S13 CONFLICT | locus intersects a `regions.by_subfam.bed` region with >= 2 names | `step3:308-339` | nested / overlapping families |
| S14 sim_ratio | single `ssearch36 -E 100` bits of the copy vs its assigned consensus / consensus self-bits | `step3:71-157` | length-free similarity; used by step 1c percentile cut (`step1c:98-116`) and plots |
| S15 arrays, twins | `tools/array_flag.py` (regular spacing vs chance null), `flankscan/fs9` flank twins (shared flanks = not independent insertions) | `SINEderella:546-547` | marks on plates, plate order; no copy removed |
| S16 boundary | step 7: extend flanks in 50 bp steps until the fraction of pairs above 45 % identity <= background + 5 points, cap 1000 bp | `step7:120-125,204-267` | sizes plate flanks; only `confirmed` is used (`step8a:130-140`) |

## Merging / clustering / extension summary

- Merging happens three times: within a window (S2, strand-aware), across
  windows per consensus (S5, strand-aware), across consensuses (S7, not
  strand-aware). Nothing clusters copies by similarity before the vote; SubFam
  chunking (`step1:247-261`) is only for proposing new consensuses from a
  30,000-copy sample and does not touch assignments.
- Extension: none at search time (no flank is kept in `extracted.fasta`);
  flanks are added only for plates (`step8a`, `boundary_refinement.tsv`) and
  flankscan windows (+-1000 bp, `flankscan/fs1_extract.sh`).
- "Candidates" are kept at three levels: `rawhits.tsv` (identity only),
  `gen-*.bed` (identity + length, post-satellite), `extracted.fasta` /
  `assignment_full.tsv` (every locus with status `assigned`, `rejected_low_bitscore`,
  `no_unanimous`). Nothing is deleted: excluded satellite hits go to
  `satellites/excluded_hits.bed`, `.before_satellites` copies are kept.

## What the thresholds mean for a replacement aligner

`notes` (my inference). A drop-in replacement for A1 must deliver, per window
or per genome region: strand, subject coordinates, % identity of the best local
alignment, and a score comparable across copies (bits), such that after
strand-aware merging a copy >= 65 % identical over >= 80 % of the consensus is
found. The vote (B1) needs a score whose ranking between consensuses is stable;
`-z 11` noise is currently absorbed by 10 repeats, so a deterministic scorer
would make the 10 cycles redundant (but `MANUAL.md:1212` presents the noise as
a feature that removes borderline copies). The 0.45 threshold is on summed raw
bits, so any score rescaling changes its meaning.
{% endraw %}
