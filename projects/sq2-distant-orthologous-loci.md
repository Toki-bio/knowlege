---
title: sq2 distant orthologous loci pilot
type: notes
---

> Type: notes

# Distant sq2 orthologous loci — pilot plan

Repo / tool context: [Toki-bio/SINE_orth_loc](https://github.com/Toki-bio/SINE_orth_loc)

This is a **project plan**, not a result. The goal is to test whether any sq2 insertion loci survive deep squamate divergence as recoverable orthologous sites, and only then ask whether they can be found in tuatara.

---

## Scientific question

Not: map all sq2 between two lizards by SINE similarity.

Rather:

> Do any **sq2 insertion events** predate deep squamate splits, such that the **same genomic locus** (unique flanks + insertion state) can still be recovered between distant taxa — and if so, can those survivor loci also be recovered in tuatara?

Keep three questions separate:

| Question | Evidence |
| --- | --- |
| Same SINE family? | sequence / consensus similarity |
| Same insertion event / locus? | unique flanks + shared or empty site |
| Same genomic neighborhood? | optional synteny support |

Only the middle one is the core claim of this pilot.

---

## Recommended pilot design

### Phase 1 species pair

**\*Gekko gecko\* × \*Heloderma charlesbogerti\***

Why this pair:

* both have high-quality chromosome-level HiFi assemblies;
* phylogenetically very distant (Gekkota vs Toxicofera / Helodermatidae);
* both have decent sq2 counts;
* if any loci survive here, they are candidates for genuinely ancient insertions.

Why not first:

* **Python / snakes** — more distant TE ecology confounders before the method is proven;
* **Rhineura** — highly rearranged burrower genomes; flank mapping may fail mechanically even if orthology exists.

### Query set

Start from **high-confidence sq2** only (e.g. sq2m / best SINEderella subset), not the full ~70k dump. Document the exact BED/source used.

### Phase logic

```text
Phase 0 (optional): calibrate on a closer species pair
Phase 1: Gekko × Heloderma — can anything survive ~220 Ma?
Phase 2: take survivor loci only → targeted search in:
         Sphenodon punctatus (tuatara)
         Dibamus, Lanthanotus, Eublepharis, Anolis, Python, ...
```

Do **not** genome-wide blast relaxed queries into tuatara first. Use Phase-1 survivors as anchors.

---

## What `SINE_orth_loc` already gives you

The existing orthology logic is already the right biology:

* flank-based locus search;
* **PM** = present in species 1, absent in species 2 (empty site in 2);
* **MP** = present in species 2, absent in species 1;
* **SINE** = present in both (shared orthologous insertion).

For ancient shared insertions, Phase 1 mainly hunts **SINE/SINE**.  
PM/MP become especially valuable later for dating within squamates.

Caveats from existing practice:

* cluster IDs are run-local, not global — prefer genomic coordinates;
* SINEderella subfamily labels and de-novo tribe labels are different systems — do not mix them when interpreting survivors.

---

## Parameter relaxation for ~220 Ma (proposal, unproven)

Claude / working discussion suggested lowering conservation cutoffs roughly like:

| Parameter idea | Recent / standard | Proposed ancient |
| --- | --- | --- |
| Flank % identity | ~65%+ | ~40%+ |
| SINE % identity | ~65%+ | ~50%+ |
| Identical flank bp floors | higher | much lower (e.g. tens of bp) |
| Flank length floors | ~150 bp | keep structural length |

This may be necessary at deep divergence. It is also the main false-positive risk.

### Required guards if thresholds are relaxed

1. **Negative control** — shuffled / scrambled target (or equivalent). Expect ~0 “orthologs.”
2. **Positive sanity** — same pipeline on a closer pair still recovers many loci.
3. **Complexity filter** — short “identical islands” must not be poly-A / SSR junk.
4. **Manual review** of the first 10–20 hits (alignments / simple plots).
5. Prefer hits with recognizable insertion architecture (SINE-sized insert; TSD if recoverable), not %ID alone.

Without the scramble control, relaxed ancient thresholds are not trustworthy.

---

## Method notes worth keeping

### Reciprocal evidence

Map both directions (Gek→Hel and Hel→Gek). Prefer unique / near-unique flank placement with consistent strand and biologically sensible presence/absence state — not merely the highest mediocre bitscore.

### Multi-scale flanks

Trying ~300 / 600 / 1200 bp can diagnose where unique mapping collapses. Use as a ladder, then freeze one primary setting for the pilot.

### Synteny

At 220 Ma, local rearrangements happen. Treat conserved neighbors as a **soft flag / ranking aid**, not an early hard kill filter.

### Masking

Pre-masking other repeats may help flanks land uniquely, but over-masking can destroy unique sequence. Treat masking as a sensitivity experiment.

### Tuatara

Ultimate stress test / outgroup to squamates. Only search **Phase-1 survivor loci** there. A hit in tuatara would be extraordinary; therefore false-positive control is mandatory.

---

## Success / failure definitions (freeze before running)

### Success (example)

≥ N loci that:

* place uniquely (or near-uniquely) by flanks in both species;
* have coherent SINE/SINE or clean empty-site logic;
* pass complexity filters;
* fail on scramble control;
* survive spot manual review.

Exact N should be chosen before looking at results.

### Failure

Scramble control yields similar counts, or no hits survive manual review / architecture checks.

Only after success: expand the survivor set to tuatara and other taxa.

---

## Non-goals for Phase 1

* full dating of all sq2;
* genome-wide Rhineura / snake mapping;
* claiming orthology from SINE % identity alone;
* relaxing thresholds until something appears, then stopping;
* mixing unrelated SINEderella label systems as if they were the same.

---

## Bottom line

```text
keep the orthology definition strict (flanks + insertion state)
relax sequence conservation cautiously
prove the relaxation with negative controls
expand phylogenetically only with survivor anchors
```

The danger at ~220 Ma is not missing a clever algorithm. It is lowering filters until non-orthologous noise looks like antiquity.
