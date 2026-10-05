---
title: TinT — transpositions in transpositions, and nested SINEs
type: notes
---

> Type: notes

# TinT (transpositions in transpositions) and SINE-in-SINE loci

Kriegs JO et al. 2007 *BMC Evol Biol* 7:190 ("Waves of genomic hitchhikers…", Galliformes CR1);
Churakov G et al. 2010 *BMC Evol Biol* 10:376 (web TinT, primate Alu chronology). Full texts and the
web tool (retrogenomics.uni-muenster.de) were not reachable when this page was written; the method below
is assembled from abstracts, search snippets of the methods, and the published model assumptions, and the
likelihood is a **reconstruction**, not a transcription of the original supplement.

## Principle

Young, active elements insert into old elements, never the reverse in time. Every nested pair
(element A inside element B) is a directional observation: A is younger than B. Counted over thousands
of nestings, the matrix *n[A][B]* (A inserted into B) orders element types (families, subfamilies) in
relative time without any molecular clock.

## Detection (RepeatMasker-based, as published)

Nested clusters are extracted from RepeatMasker coordinates with five criteria; three are known verbatim:

1. fragmented and inserted elements lie on the same query sequence;
2. *stringent*: both parts of the fragmented (host) element share the RepeatMasker element index and the
   nested element's index is higher; *relaxed*: the same class for the host parts is enough;
3. every element (start/end) is ≥ 20 nt.

The other two concern adjacency/orientation (not recovered). Output: an insert × host frequency matrix
plus copy numbers per type.

## Model (TinT assumptions → closed form)

Published assumptions: each type has **one** period of activity, **normally distributed**, and **no target
site preference**. Then for types A, B with activity times *T_A ~ N(μ_A, σ_A²)*, *T_B ~ N(μ_B, σ_B²)*, an
insertion of A can land in B only if B already exists, so among nestings between A and B:

P(A is the insert) = P(T_A > T_B) = Φ((μ_A − μ_B) / √(σ_A² + σ_B²))

— a Thurstone paired-comparison model. With a common σ (case V) the μ are identifiable up to a shift and
are fitted by maximum likelihood over all off-diagonal counts; bootstrap gives intervals. Display as in
TinT: one bar per type, centre = μ, ellipse ≈ 75%, line ≈ 95% of the activity distribution.

Caveats that follow from the model:

* pairs observed only one way (A always inside B) give a **bound**, not a value; report it as such;
* σ per type is weakly identifiable from nesting counts alone (self-nesting n[A][A] carries some width
  information); case V is the honest default;
* target preferences (e.g. Alu into Alu A-rich tails), deletions of old hosts and annotation bias towards
  young, intact elements all distort counts.

## Implementation for SINE_orth_loc (`sine_nest.py`)

* `scan`: nhmmer (HMMER) instead of RepeatMasker; pieces ≥ 20 nt; a nesting needs two host pieces whose
  **consensus coordinates join up** (overlap = TSD) around the insert, and the **target site duplication is
  verified in the sequence** (direct repeat 6–25 bp flanking the insert); host pieces that nhmmer extends
  into the (similar) insert are trimmed back first — without this, 26% of simulated nestings were missed.
* Also classified: split annotation (no insert), dimers, satellites (periodic units with shared spacers —
  not independent insertions), close copies.
* `tint`: the Thurstone fit above; on simulated insertion histories (6 types, uniform targets) the true order
  is recovered with Spearman 1.00, also with unequal activity widths and ~180 events.

Results: simulated genomes — 106/110 nested inserts found, 0 false, TSD in all, host the more diverged
element in all. Real *D. valentini* / *D. mixta* HiFi assemblies (dar_squam1): 47 / 32 nestings; host more
diverged in 41/47 and 30/32 (median identity to consensus 81–84% vs 93–94% for inserts), TSD found in
36/47 and 25/32 (median 14 bp). Same-family SINE-in-SINE is rare there; subfamily-level TinT needs the
SINEderella subfamily consensi as the library.

## Why nesting matters for orthology

The young insert's flanks *are* the old host, and in a window holding both, the consensus aligns to the
younger, more similar element: a flank-based pipeline can pair an old host in one genome with a young
insert in another. In simulation this produced 2–6 wrong host and 2–4 wrong insert calls per genome pair;
comparing the compound locus through its outer flanks and testing each element by junction vs empty-site
junction removed them (0 wrong) and raised recall of nested inserts from ~50% to ~95%. A nested pair is
also a phylogenetic character with a built-in age order: the insert is younger than the host in every
lineage.

Related: [TE-insertion anchor graph](https://toki-bio.github.io/knowlege/projects/te-anchor-graph-pansineome/),
[SINE_orth_loc](https://github.com/Toki-bio/SINE_orth_loc).
