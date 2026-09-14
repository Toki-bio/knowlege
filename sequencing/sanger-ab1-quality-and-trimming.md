---
title: Sanger ab1 quality assessment and trimming, done right
type: explanation
---

> Type: explanation

# Sanger ab1 quality assessment and trimming, done right

## Short answer

**Never trust the basecaller's Phred score or the plain letter sequence alone**, and never trust a rendered chromatogram plot glanced at by eye either. Both fail in ways that look confident. The only defensible approach is: extract the raw four-channel trace numerically, measure peak topology per base, and cross-check the result against a real reference (BLAST) before deciding what's usable gene sequence vs. noise, primer carryover, or an unresolved blob.

This note exists because both failure modes were hit directly, in the same afternoon, on the same file, working through the [Artemia_parthenogenetica](https://github.com/Toki-bio) COI/18S/ITS resequencing.

---

## Failure mode 1: Phred score can be confidently wrong

Phred/basecaller quality is derived from the *relative* shape and separation of the four channel traces at a position, not their absolute height. When a read runs out of real signal (weak template, failed extension, run-off past the true product), all four channels can sit near flat baseline noise while still showing enough relative wiggle for the caller to keep assigning respectable Phred scores for hundreds of bases past the point where nothing real is left.

Observed directly: two ITS reads (`G4`, `H4`) had raw trace amplitude collapse to <150 RFU after roughly 2,000 of ~14,000 total scan positions — yet Phred score in that dead zone (~26 average) was *higher* than in the genuinely good signal region earlier in the same read (~19-23 average). A plain "trim below Q20" pass keeps the garbage and cuts into the real data. Quality-based trimming (Mott's modified algorithm, or a sliding-window cutoff — see below) is still useful as a first pass, but it is not sufficient on its own.

## Failure mode 2: eyeballing a rendered trace plot is also unreliable

The natural next step — plot the four channels and look at it — is necessary but not sufficient either. Two concrete misses from the same 287bp 18S read (`C1`):

* A region initially read as "recovers to clean signal a few bases later" was, on actual per-base comparison to a BLAST-matched reference, still wrong at 5 of 10 positions.
* A region at the read's 3′ end that visually looked like a total signal collapse ("noise, discard") turned out, once BLASTed, to match the real 18S reference at >90% identity with only two or three small indels. It was real, low-amplitude-but-genuine tail signal, not noise.

Both errors happened because a rendered plot was described qualitatively ("looks messy", "recovers") instead of measured. Rendered images compress and smooth exactly the distinction that matters: a genuinely resolved single peak vs. a wide, flat-topped merged peak (a "blob") that the basecaller sliced into several fake single-base calls, each of which can still print a confident, unambiguous letter.

---

## What actually works: numeric peak-ratio method + BLAST cross-check

### 1. Per-position peak measurement (not a visual read)

Using Biopython's `SeqIO.read(path, "abi")`, the four raw channel traces live in `annotations["abif_raw"]["DATA9"..."DATA12"]` (channel-to-base order from `FWO_1`), and per-base scan positions in `PLOC1`/`PLOC2`.

For every called position `i`:

1. Take the trace window between the midpoints to its neighboring calls — `left = (ploc[i-1]+ploc[i])//2`, `right = (ploc[i]+ploc[i+1])//2` — this is that call's actual slot in the trace, not an arbitrary fixed window.
2. `heights = {base: window[base].max() for base in 'ACGT'}`.
3. Rank channels by height; `ratio = second_highest / highest`.
4. Check whether the dominant channel has its **own local maximum** inside this slot (e.g. `scipy.signal.find_peaks`), or whether it has none — the fingerprint of a slot that's really just an arbitrary slice of a neighboring blob peak.
5. Flag the position if `ratio ≥ threshold` OR there's no local apex of its own OR the top height is near that file's own noise floor.

**Threshold — this is an established parameter, not something to improvise.** `sangerseqR::makeBaseCalls` (the standard R tool for exactly this problem) uses **0.33** as its default secondary/primary cutoff for calling a position heterozygous/ambiguous. Forensic mtDNA heteroplasmy analysis uses a tighter **10-20%**, because there a missed real minor variant is costlier than a false positive. Use 0.33 as a sensible default; tighten to 0.2 if the downstream use can't tolerate missing a real secondary peak.

Collapse flagged positions into contiguous runs. Isolated single-base flags scattered through an otherwise clean stretch are ordinary Sanger noise-floor wobble, not a problem region — what matters is multi-base contiguous runs. A position can fail this test even when the basecaller printed a confident, unambiguous single letter; that mismatch between basecaller confidence and actual peak topology is exactly the thing worth catching.

### 2. BLAST cross-check, not trust in the trace alone

Independently BLAST the raw, full-length called sequence against NCBI `nt` (the URL API — `CMD=Put` then poll `CMD=Get` — works fine from a plain `curl` script, ~15-30s turnaround). Diff the returned query/subject alignment position by position, not just the headline percent identity.

A stretch is trustworthy gene sequence only when **both** signals agree: real peak topology (method above) *and* it falls inside a BLAST alignment matching a real reference at high identity. This distinguishes two categories that look similar at a glance but aren't:

* A stretch **outside** the BLAST alignment (commonly the first few bases of a read) isn't "bad quality data" — it's not the target locus at all, typically primer/adapter carryover before the polymerase is reading real template.
* A stretch **inside** the alignment with clean peak topology but only ~50% identity to every top hit is grounds to suspect contamination or the wrong locus, not a trimming problem.

### 3. Coarse algorithmic trim as one input, never the verdict

Two standard trimming approaches, both implemented in `sangeranalyseR`:

* **M1 — modified Mott algorithm** (also in Phred, Biopython): cutoff on a log scale — Q10 → 0.1, Q20 → 0.01, Q30 → 0.001. Trims contiguous low-quality runs from both ends.
* **M2 — sliding window** (Trimmomatic-style): default window size 10, cutoff Q20.

Run one of these for a coarse first pass, then verify the result against the numeric peak method and the BLAST alignment. Treat its output as a proposal.

---

## Practical checklist

1. Extract raw traces + Phred + base calls (Biopython `abi` parser).
2. Plot the four channels + a Phred bar underneath, at a width where individual peaks are actually legible — this is for orientation, not for issuing a verdict.
3. Run the numeric peak-ratio method on every position (threshold 0.33 default).
4. BLAST the raw called sequence against `nt`; diff the alignment position by position.
5. Run a coarse quality trim (M1 or M2) as one more input.
6. Report contiguous flagged runs with the actual numbers (ratio, own-apex yes/no, BLAST match/mismatch), and state a real opinion on which stretches are gene signal vs. artifact — backed by the evidence, not hedged into "you decide" when the evidence already supports a conclusion.
7. Only then write the final FASTA — lowercase or IUPAC-coded for low-confidence-but-called bases, discarding what neither method supports.

## The ratio test can itself be wrong on a globally weak read

The ratio test assumes background noise is small relative to a real peak. On a **globally low-amplitude read** — weak template or dilute reaction, where the whole trace sits at low absolute RFU throughout (real peaks of 40-150 RFU against ~20-50 RFU background), rather than a read that is dying partway through — that assumption breaks. Ordinary baseline wobble becomes a large fraction of every peak height, and the ratio test flagged roughly half of two ITS reads in this project as unreliable, flatly contradicted by BLAST: those same reads were ~95-97% identical to the matched reference across nearly their full length (only ~24 real mismatches across ~800 aligned bases). **Before converting a flagged run to N, check the actual local BLAST mismatch/gap density in that run.** If the reference disagrees at only a handful of positions in a "flagged" stretch, the numeric test was picking up this file's noise floor, not a real failure — downgrade to lowercase instead of N. Reserve N for runs where BLAST-confirmed disagreement density is genuinely high, or where the run has no BLAST coverage at all. The BLAST cross-check is load-bearing for the N-fill decision itself, not just for locating the read's overall start and end.

## What to watch for, by artifact type

* **Primer dimer / injection artifact** — noisy, low-amplitude region at the very start, before the polymerase is reading real template. Falls outside the BLAST alignment.
* **Dye blob** — broad, flat-topped peak spanning several nominal base-call slots, with no real local maximum for most of those slots. Flagged by "no own apex" in the numeric method above, regardless of how confidently the basecaller lettered it.
* **Spectral pull-up / oversaturation** — amplitude above roughly 10,000-15,000 RFU bleeds into neighboring channels, mimicking a double peak. Check whether the "secondary" peak is just bleed-through from a huge immediate neighbor before calling it ambiguous.
* **Real 3′ signal decay** — amplitude drops toward the end of a read, but if it still resolves single dominant peaks per slot and BLASTs cleanly to the reference, it's real, just weak — not noise to discard by reflex.

## Bottom line

None of this is novel. Phred's own 1998 definition is built from peak spacing, resolution, and uncalled/called peak ratio — this is a simplified, transparent reconstruction of the same category of signal, made explicit because the finished Phred score alone hides exactly the distinction that matters. The failure isn't a missing method; it's substituting a plausible-sounding description of a plot for the numbers that were already sitting in the file.
