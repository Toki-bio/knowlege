---
title: Composite vs unit SINE assignment
type: explanation
---

> Type: explanation (measured results marked; one design choice marked as notes)

# Composite vs unit: assigning a SINE copy when the bank has both a dimer and its monomers

## The problem

A consensus bank can hold a composite element (a dimer, or unit + linker + unit, such as Alu vs FLAM/FRAM, B1 vs
PB1/B1-dID, or the rsi bat composites r1_r3_P1 = r1 + 39 bp + r3) next to its units. A genomic copy can then be:

- a **single unit** (free monomer),
- a **full composite**,
- a **truncated composite** (often 5'-truncated, so it looks like the downstream unit, maybe with a bit of linker).

For a single-unit copy, a local aligner (ssearch36, BLAST, cross_match) aligns the **same bases at the same identity**
to the unit and to the unit's half of the composite. The scores are equal in principle, so "best score wins" cannot
decide.

## What existing tools do (survey, 2026-10-16)

No published tool decides unit vs full composite vs truncated composite per genomic copy. Each one keeps the
best-scoring consensus, and joins fragments only within one name or a declared name pair. Report:
`Raks_COI/reports/Composite SINE annotation tools.md`.

| tool | relevant logic | borrowable |
|---|---|---|
| RepeatMasker (Perl, OSL-2.1) | Alu vs FLAM/FRAM settled by **refinement**: a re-alignment counts only if it covers the copy to within 10 bp, so a full Alu is never renamed FLAM. Only truncated-composite rule found anywhere: `ID_B1` hit ending before consensus position 90 is renamed `ID`. Fragment joining needs the same or substring-related names. | idea (reimplement, copyleft) |
| LTR_retriever `solo_finder.pl` (Perl, GPL-3) | keeps the full LTR-int-LTR consensus **plus a table of where each LTR sits in it**; a hit is a solo LTR when it stays inside the LTR interval (covers >= 0.8 of it, spans <= 1.2x) and no internal region is near | yes: the coordinate rule |
| dfamscan.pl | drops a lower hit when >= 75 % of it lies inside a higher one; no composite logic | no |
| RepeatCraft, OneCodeToFindThemAll, REPET long join, REANNOTATE, TEnest | fragment joining by name / name dictionary (LTR + int); no unit-vs-composite decision | partly (name dictionary = declared composition) |
| MCHelper | flags a consensus "chimeric" when a terminal direct repeat is > 35 % of it; would flag real dimers | as a check only |
| AnnoSINE_v2, SINE_Scan, SINEBase | no composite assignment; SINEBase describes families as modules (head, body, tail) for one query at a time | no |

Field practice in SINE papers: count a copy as composite only at >= 90 % coverage of the composite consensus
(Kosushkin et al. 2024, Urop; Ray et al. 2015, Ves), and accept a free monomer as independent only with its own TSD.

## What happens in a ssearch36 voting pipeline (SINEderella, measured)

SINEderella step 2 votes 10 times with `ssearch36 -z 11` (shuffle statistics). Measured on a toy (two families, unit A
120 bp, unit B 150 bp, dimer AB = A + 25 bp + B; 40 copies per class, 5-15 % diverged; 2 Mb genome):

- **The tie is noisy, not exact.** `-z 11` estimates statistics per consensus, so the same alignment (87.65 % over
  162 bp) scored 149-159 bits against the unit and 137-150 against the dimer, varying by cycle.
- **Result on main (296fa33):** dimers were always called dimers, but single units flipped run to run: one family's
  A units were 99/100 firm in one run and 0/40 in the next; the B units, identical to the dimer's B half, were never
  firm (no_unanimous). Nothing was mislabelled; copies were lost.
- **Truncated dimers often lose their head in step 1**, before any assignment: the dimer's hit is shorter than
  0.8 x the dimer length, so sear drops it, and the locus is only the B part.

## Fix that worked: declared composition + coordinate rule

1. A table declares each composite and the **intervals of its units on the composite** (`r1_r3_P1  r1_9seqs:6-161,
   r3_58seqs:201-401`). Intervals come from the composite's build layout (flankscan peaks.tsv), placed in order
   by aligning the declared unit parts (`tools/place_modules.py`).
2. In each voting cycle, only for a declared unit u and its composite c within 15 % in bitscore: the copy's alignment
   to c decides. If it reaches >= 10 bp beyond every interval of u into another module or the linker (on a side where
   c has > 30 bp of other sequence), it is c. Otherwise it is c only if > 2 identity points higher, else u.
3. No table: behaviour identical to before.

Toy result: 316/320 correct in each of two replicates, 0 units called dimer, 0 dimers called unit.

## Two things that failed (keep them out)

- **"Explains more of the copy" without module coordinates** moved 643 single r6 copies to the composite r5_r6_P26:
  P26's built-in copy of r6 has a ~15 bp longer 3' tail, so a tail difference looked like composite evidence.
- **Detecting unit-inside-composite pairs from the bank automatically** paired related length versions and
  composites too (r9 inside r7, P2 inside P1) and moved 448 firm copies between ordinary subfamilies. Units that
  share a head (tRNA-derived SINEs) also align into the wrong module: r8 onto P18's r10 head, r5 onto P26's r6 part.
  Composition must be **declared**, and intervals placed **in declared order**.

## rsi (Rhinolophus sinicus) finding

Bank of 21 (10 peel units r1-r10, 4 MEG, 7 composites), 67,119 copies, step 2 only, same step-1 copies
(therioserver `/home/toki/rhin/rsi_tiebreak_20261016/`). Noise baseline (main vs main): 683 firm <-> not firm, 0 firm
moves. Declared composites with explicit intervals (branch `coverage-tiebreak` d7c09b9): 0 firm moves between ordinary
units/subfamilies; 421 composite-related firm moves.

- **r1_r3_P1 (14,500 copies) is not absorbing single r3 or r1.** With correct intervals, P1 changes by +1, r3 by +13,
  r1 by +3. The P1 copies carry the r1 + linker + r3 junction; the earlier suspicion of absorption was wrong for P1.
- **r6 -> r5_r6_P26 (354 firm copies):** 238 of them align >= 10 bp (median 33) into P26's linker / r5 part upstream
  of r6, so they look like 5'-truncated P26 copies that main called r6 (r6's bitscore was ~10 % higher because the
  extra part aligns at lower identity). Needs judgement on the alignment.
- **r10 vs r10_r8_P18:** P18's second element is not r8 (r8 aligns only to P18's shared tRNA head, positions 2-133);
  P18 = r10 + 276 bp of other sequence. r10's real 3' tail runs ~21 bp past the r10 consensus (stage 8 EXT3: r10 +21),
  and P18 contains that tail after its r10 part, so plain r10 copies "reach" into P18: 264 copies vote P18 10/10 and
  fail P18's bitscore threshold. A short unit consensus makes its own tail look like composite evidence: the declared
  interval must include the unit's real tail (here r10:6-171).
- Module placement: whole units that share a head align to the wrong module (r8 onto r10's head in P18, r5 onto r6's
  part in P26, one r5 in the homodimer P48, overlapping r1/r3 in P2). Placing the declared unit parts from flankscan's
  layout (peaks.tsv: upstream unit end, linker length, downstream start) in order fixes P1, P2, P26, P34, P48.

Review alignments (ViewAlign): Tal `rsi_composites_test/` (r6 vs P26, r5 vs P48, r10 vs P18, r3 vs P1).

## Notes (not established)

- glsearch36 (global in the consensus) would penalise the unaligned part of a composite directly; untested.
- A free monomer could additionally require its own TSD, as in the field practice; not implemented.
