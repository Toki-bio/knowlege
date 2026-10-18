---
name: trna-similarity-scan
description: How to test whether a SINE consensus (or a part of a composite SINE) has a tRNA-, 5S rRNA- or 7SL-derived head, and from which tRNA isotype, honestly - including degraded heads that tRNAscan-SE and Rfam miss. Covers pairwise search against tRNA gene sets with a dinucleotide-shuffle null, SINEBase's 60 %/60 nt rule, Infernal/Rfam and tRNAscan-SE settings and thresholds, Pol III box scans, folding/covariation as support only, pitfalls, and a tested recipe with rsi (Rhinolophus sinicus) results. Use when classifying SINE heads, checking composite SINE parts, or reviewing any "tRNA-related" claim.
---

# tRNA / 5S / 7SL similarity of SINE heads

> Written 2026-10-18 from a deep-research report (`references/report.md`, sourced) and one measured run (rsi).
> Knowledge page: knowlege `methods/trna-similarity-scans-for-sines.md`.

## Rules (from the sources)
- The field's test is **similarity to a real RNA gene**, not motifs: SINEBase (Vassetzky & Kramerov 2013) calls a head
  tRNA-, 5S- or 7SL-derived at **>= 60 % identity over >= 60 nt**; the isotype is the best-matching tRNA.
- **A / B Pol III boxes are "an indication, not a proof"** - 5-11 nt patterns match by chance; 5S-derived heads have
  no B box (A box + IE + C box), so a B-box requirement misses them (megabat MEG-R* families).
- **Default detector thresholds are built to keep SINEs OUT of tRNA gene lists**: tRNAscan-SE reports at 20 bits, its
  high-confidence filter needs isotype score 70; Rfam gathering thresholds tRNA 29, 5S 38, 7SL 70 bits (a 7SL-derived
  head covers part of the 299-nt model and falls below 70 - rank by E-value instead).
- **Anticodon calls of decayed heads are noise**; the isotype-specific model (tRNAscan `--detail`, "IPD" = anticodon
  and isotype disagree) and the best-hit isotype are the signal.
- **Every score needs a null**: dinucleotide-preserving shuffles (`esl-shuffle -d`, `ushuffle -k 2`) and, ideally,
  random genome windows of the same length and GC. Report "not detected at threshold", never "not tRNA-derived".
- Folding (RNAfold, RNAalifold) and covariation (R-scape) only support a call; old heads keep the boxes, lose the pairing.

## Tested recipe (rsi, 2026-10-18; `references/trna_rsi.sh`, `references/summary.py`)
1. Queries: the consensuses plus every stretch of each composite (SINEderella `composite_units.tsv`), so each unit's
   head is tested on its own (a single best hit per consensus hides the second head of a dimer).
2. `ssearch36 -3 -m 8 -E 10 queries.fa tRNAs.fa` against host-near tRNA genes (human GtRNAdb hg19, 419 genes - the hg38
   download was truncated, check counts; Myotis lucifugus is the only bat in GtRNAdb). Keep best tRNA, identity, length,
   query span, bits; best hit of a different amino acid = margin.
3. Null: `esl-shuffle -d -N 20 --seed 42 queries.fa` -> same search; a head counts when its bits exceed the best shuffle.
4. `cmscan -E 1 --tblout cm.tbl Rfam.cm queries.fa` (all Rfam: tRNA, 5S, 7SL, ...).
5. `tRNAscan-SE -E -H --detail -X 10 -o out -f ss queries.fa` (isotype scores; pseudo / IPD flags).
6. Tabulate (summary.py): SINEbase 60/60 pass, bits vs null, 2nd isotype, Rfam, tRNAscan, A/B box positions.

## Measured on rsi (therioserver `/home/toki/rhin/rsi_main_20261018/trna/`, Tal `rsi_main_20261018/trna/`)
- Every r-unit: tRNA-Ile head, 73-93 % over 66-75 nt, bits 39-82 vs shuffle max <= 29; next amino acid 8-37 bits lower.
- **Rfam missed r2, r3, r4, r6; tRNAscan-SE missed r3, r4**; step 2 + 3 found all (the 2026-09-29 run needed Dfam's
  online tRNA-Ile HMM for them).
- Composites: one head per unit (P1 two, C11 three, P26/P34/P48 two, P18 one); linkers and tails: none.
- Controls: MEG-RL/RS 5S (Rfam 1-119, E 1e-24), MEG-T2 tRNA-Val, MEG-TR tRNA-Val + 5S.

## Locations
- therioserver `/home/toki/refs/smallrna/`: `Rfam.cm` (+ index), `hg19-tRNAs.fa` (419, complete), `hg38-tRNAs.fa`
  (TRUNCATED, 123 - do not use as a full set); conda env `rnatools` (cmscan, tRNAscan-SE), `sinederella` (ssearch36, esl).
- Dfam curated consensi: therioserver `/home/toki/refs/dfam_current/`.
