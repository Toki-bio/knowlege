---
title: Alu graph-map SV proposal evaluation
type: notes
---

> Type: notes

# Evaluation: Alu-anchor graph for cheap cancer SV profiling

This is a critical evaluation of a proposal to treat Alu loci as **fiducial markers** on a probabilistic graph and measure cancer structural variation (SV) by changes in Alu presence / inter-Alu spacer topology — using short reads, outward Alu-PCR, and ONT adaptive sampling.

Status: **interesting research / translational concept**, not a validated assay. Strengths and landmines are separated on purpose.

Related in this knowledge base:

* [Alu-anchor graph topological profiling](https://toki-bio.github.io/knowlege/projects/alu-anchor-graph-topological-profiling/) — the working **concept spec** (graph, \(W_u\), \(H_\text{norm}\), COLO829 PoC). This page is the **evaluation**.
* [Alu ddPCR somatic burden screen](https://toki-bio.github.io/knowlege/projects/alu-ddpcr-somatic-burden-screen/) — bulk young-Alu **activity**, not topology.
* [TE-insertion anchor graph — pan-SINEome and Alu topology](https://toki-bio.github.io/knowlege/projects/te-anchor-graph-pansineome/) — the same graph as a general extension of the pan-SINEome; further critique (copy number, low-pass CNV baseline, matched normal).

---

## The proposal in brief

### Graph model

Model the genome as a probabilistic variation graph \(G = (V, E, P)\):

* **Vertices \(V\)** — Alu loci (not treated as identical). Annotate:
  * coordinates (GRCh38 + pangenome alt loci),
  * subfamily / divergence profile (AluJ / AluS / AluY …),
  * uniqueness weight \(W_u \in [0,1]\) from k-mer uniqueness  
    (old divergent AluJ → high \(W_u\); young AluY → low \(W_u\)).
* **Edges \(E\)** — inter-Alu spacers with sequence + mappability.
* **Probabilities \(P\)** — trained on a ~100-genome pangenome:
  * \(P(v_i\ \mathrm{present})\),
  * spacer-length distributions (mixture models preferred over a single Gaussian).

### Technology stack

| Layer | Role | Strength | Blind spot |
| --- | --- | --- | --- |
| Multiplex short reads on “Alu-tail + spacer” index | junction / presence-absence for high-\(W_u\) nodes | cheap, scalable | weak long-range topology |
| Outward Alu-PCR | enrich edges (spacers) | makes spacers visible to short reads | length bias; misses long edges |
| ONT low-pass (1–3×) | span \(v_i \to e_{ij} \to v_j\) | proves adjacency + length | cost / throughput if genome-wide |
| ONT Read Until adaptive sampling | enrich Alu/spacer-bearing molecules | turns cheap run into high on-target depth | reference design is critical; poor for ultra-short cfDNA |

### Dual clinical outputs

1. **Global topological entropy / “CIN score”** — genome-wide deviation of patient graph from expected normal graph → agnostic chromosomal-instability biomarker.
2. **Hotspot subgraph tests** — Bayesian / CI tests on driver neighborhoods (BRCA1, EGFR, MYC, …) → actionable calls.

Pitch paragraph theme: reduce SV detection from “deep WGS of 3 Gb” to “measure distances/adjacencies among ~10⁶ Alu landmarks.”

---

## What is genuinely strong

### 1. The conceptual reframe is real

Treating Alus as **landmarks rather than noise** is a legitimate sparse-coordinate idea. The genome already has ~1 million Alu insertions; old divergent copies behave almost like unique sequence. Weaponizing that abundance for topology is the proposal’s core insight.

### 2. \(W_u\) is not cosmetic — it is load-bearing

Without a divergence / uniqueness weight, short-read mapping collapses. Distinguishing:

```text
AluJ / high divergence  → usable anchor
AluY / near-identical   → ambiguous multi-mapper
```

is exactly the right first-order correction. This is where Alu biology expertise matters.

### 3. Mixture models for spacer lengths are the right default

Healthy genomes already contain common indels and mobile-element variation inside spacers. A single Gaussian will either:

* overstate variance and kill sensitivity, or
* understate multimodal haplotypes and create false SV calls.

Mixture / haplotype-aware edge models are appropriate.

### 4. Dual output (entropy + hotspots) is strategically smart

CIN-like global scores and actionable local calls serve different customers. That dual framing is fundable and scientifically coherent — **if** entropy can be de-noised.

### 5. Adaptive sampling is the right hardware bet for solid tumors

If Read Until can keep molecules that start in Alu / near-Alu sequence and reject deep unique sequence without Alu context, on-target depth can rise dramatically without deep WGS cost. That is the commercial hinge.

---

## Critical landmines (take these seriously)

### 1. Young / polymorphic AluY nodes will break naive traversal

Nodes with low \(W_u\) and \(P(\mathrm{present}) \approx 0.4\) create branching ambiguity. Short reads alone will not resolve them. Design rule:

```text
high-Wu nodes  → short-read / PCR primary evidence
low-Wu nodes   → ONT-required, or drop from short-read graph
```

Do not pretend all Alus are equal anchors.

### 2. Adaptive sampling reference catch-22

If the reject/keep reference includes **full wild-type spacers**, a translocation that keeps an Alu but switches the spacer to another chromosome may be rejected as “not matching.”

**Necessary fix:** adaptive-sampling targets should emphasize:

* conserved Alu subfamily consensuses / tails,
* short immediate flanks,

**not** entire expected wild-type spacer sequences.

Keep-if-Alu-like; do not require the healthy spacer to continue.

### 3. Outward Alu-PCR is deletion-biased

PCR competition favors short products. Expected behavior:

* good for shortening edges (deletions),
* bad for long insertions / large duplications on that edge,
* also prone to Alu–Alu artifactual chimeras.

Treat Alu-PCR as an enricher for a subset of edges, not as topology ground truth.

### 4. Entropy vs technical chimera / library artifact

High “graph entropy” can come from:

* true chromothripsis / BFB / CIN,
* PCR chimeras,
* ligation artifacts,
* mapping errors at low-\(W_u\) Alus,
* tumor purity / subclonality.

Match-normal baselines are mandatory. Without them, CIN scores will not survive review.

### 5. Liquid biopsy claim is currently oversold

cfDNA is ~160 bp. Read Until needs enough early bases to decide reject/keep; ultra-short fragments are largely sequenced before rejection helps. Safer roadmap:

```text
v1: solid tumors (+ maybe SV-driven hematologic malignancies)
later: cfDNA using short-read junctions / dPCR, not ONT adaptive sampling
```

### 6. Alu density is uneven — coverage of SV space is uneven

Alus are not a uniform lattice. Some clinically important regions are Alu-rich; others are Alu-poor. Your sensitivity will be **locus-dependent**. Hotspot design must start from Alu density maps, not only from gene names.

### 7. Scale and training reality

~10⁶ nodes × edges is tractable but not trivial. “100-genome pangenome” baselines should start from public graphs (**HPRC** / related) rather than inventing a private 100-genome stack first. Still expect:

* incomplete Alu annotation concordance across assemblies,
* coordinate liftover pain,
* polymorphic presence/absence as first-class state, not noise.

### 8. What this does **not** replace

It is not:

* a full WGS replacement for all mutation classes,
* a base-perfect assembler,
* a pure Alu-retrotransposition assay (that is the ddPCR burden idea),
* automatically liquid-biopsy ready.

It is a **targeted topological sensor** for a subset of SVs that perturb Alu–spacer adjacency/length.

---

## Verdict on the 3-paragraph pitch

Professionally strong. Clear problem → method → dual utility.

Suggested punch-ups:

1. Quantify the reduction: “3 Gb assembly problem → discrete graph of ~10⁶ Alu landmarks.”
2. Clarify Read Until: low-pass **input** becomes high-depth **on-target** sequencing.
3. Soften / stage liquid biopsy: lead with solid tumors; mention cfDNA only as a later short-read/dPCR path.

Overall pitch grade: **fundable concept if the first in-silico PCAWG/HPRC validation is honest about sensitivity by SV class and Alu density.**

---

## How to handle initial in-silico training (HPRC question)

Yes — start with public pangenome resources rather than waiting for a custom 100-genome cohort.

Recommended first stack:

```text
1. HPRC / public pangenome graphs + GRCh38/T2T Alu annotations
2. Build node set = Alu loci with divergence + Wu
3. Build edges = intervening unique sequence / distances
4. Estimate P(present) and spacer-length mixtures across haplotypes
5. Simulate deletions / duplications / inversions / translocations
6. Ask: which SV classes move which edges/nodes outside CI?
7. Benchmark against PCAWG / other gold-standard SV truth sets
8. Only then design wet-lab primer / adaptive-sampling libraries
```

Success metric should not be “detects all SVs.” It should be:

> For which SV classes, sizes, and genomic neighborhoods does an Alu-graph sensor recover truth cheaper than deep WGS — and where does it fail?

---

## Relation to the Alu ddPCR burden idea

Keep these as **two products**, not one:

| Idea | Question |
| --- | --- |
| ddPCR young-Alu burden | Did Alu retrotransposition burst? |
| Alu graph-map SV | Did chromosome topology scramble relative to Alu landmarks? |

A tumor can have:

* high CIN / many SVs with little new Alu activity, or
* Alu bursts with modest large-scale topology change, or
* both.

Do not market them as the same assay.

---

## Go / no-go checklist before calling it brilliant in public

```text
[ ] Wu separates usable vs unusable anchors on real pangenome data
[ ] Adaptive-sampling keep/reject library works for novel Alu–Alu / Alu–foreign junctions
[ ] Entropy separates match-normals from high-CIN tumors without chimera inflation
[ ] Hotspot sensitivity estimated vs PCAWG by SV class and Alu density
[ ] Liquid biopsy claims deferred until fragment-length reality is addressed
[ ] Cost model includes failed runs / enrichment bias / compute, not only flow-cell fantasy
```

If those pass, the proposal is more than elegant — it is developable.  
If they fail, the insight about Alus as fiducials may still survive in a narrower solid-tumor ONT product.

---

## Bottom-line evaluation

**Brilliant as a framing:** yes — Alus as probabilistic landmarks + tiered cheap topology sensing is a sharp idea.

**Brilliant as an immediately deployable clinical platform:** not yet — young Alu ambiguity, Read Until reference design, PCR bias, entropy artifacts, uneven Alu density, and cfDNA physics are real failure modes.

**Best next move:** HPRC/T2T in-silico graph construction + PCAWG SV recovery curves, with explicit stratification by Alu density and SV class. Do not lead with liquid biopsy.

That is the honest grade: **high-upside concept, validation-gated, not magic.**
