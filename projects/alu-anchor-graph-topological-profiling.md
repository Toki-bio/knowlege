---
title: Alu-anchor graph topological profiling
type: notes
---

> Type: notes

# Topological genomic profiling via an Alu-anchor graph

This is a **working concept spec**, not a result. It is the same family of idea as the [critical evaluation](https://toki-bio.github.io/knowlege/projects/alu-graph-map-sv-proposal/) — Alus as a coordinate system for cancer structural variation (SV) — written out as a graph, two scores, a COLO829 proof-of-concept, and a claimed novelty split versus an older Alu-cancer patent.

Related, different question: [Alu ddPCR somatic burden screen](https://toki-bio.github.io/knowlege/projects/alu-ddpcr-somatic-burden-screen/) counts young-Alu retrotransposition. This spec measures **whether the map between Alu landmarks has been scrambled**.

Shared data model with the SINE work: [TE-insertion anchor graph — pan-SINEome and Alu topology](https://toki-bio.github.io/knowlege/projects/te-anchor-graph-pansineome/) (junction-based anchor uniqueness, read-validated node calls, edge layer implemented on Darevskia).

Status: pre-validation. Numbers below (anchor fractions, *H*_norm fold-changes, COLO829 4–8×) are **design targets**, not measurements.

---

## The idea

The human genome already has ~1.1 million Alu insertions, dense in gene-rich DNA (order one per 3 kb). Instead of assembling a patient’s 3 Gb to find SVs, treat those insertions as **nodes** and the unique sequence between neighboring insertions as **edges**. Deletions, duplications, inversions, and translocations change distances and adjacencies among those landmarks.

Compare a sample’s observed inter-Alu spacer lengths (and missing adjacencies) to a **pangenome baseline**. Two intended outputs from cheap sequencing:

1. a genome-wide **CIN-like score** (*H*_norm) — “how scrambled is this genome relative to the baseline?”
2. **hotspot calls** on driver neighborhoods (BRCA1, EGFR, MYC, …)

The reduction, if it works, is operational: measure a discrete Alu graph, not assemble the whole genome.

---

## Graph \(G = (V, E, P)\)

### Nodes — Alu loci, not interchangeable

Each node carries coordinates (GRCh38 + pangenome alts), subfamily, divergence from consensus, strand, and a uniqueness weight \(W_u \in [0,1]\). Optional later layer: ONT native CpG methylation at the Alu.

The spec’s intended uniqueness gradient:

| Class | Rough divergence | Role |
| --- | --- | --- |
| AluJ (old) | ~25–50% | high \(W_u\); candidate **unconditional** short-read anchors |
| AluS | ~10–25% | mixed; **conditional** anchors |
| AluY (young) | ~0–10% | low \(W_u\); ambiguous; long-read or drop from the short-read graph |

Polymorphic insertions get a presence probability \(P(v_i\ \mathrm{present})\) from haplotypes, not a forced “always there.”

### Edges — spacers between adjacent Alus

An edge is the interval between two neighboring Alus (the spec caps gaps at ~50 kb for the first graph). Annotations: spacer sequence / mappability, and **orientation pair** (DD, DI, ID, II). Inverted pairs (DI/ID) are flagged as deletion-prone Alu–Alu NAHR candidates, not as the whole method.

### Baseline \(P\) — pangenome, not a matched personal reference

Trained once on ~100 HPRC haplotypes:

* node presence frequencies;
* edge length as a distribution (the spec writes a Gaussian \(\mu, \sigma\) with a 3σ flag and a 5 bp floor on \(\sigma\)).

Healthy common indels are supposed to sit inside that interval; only outliers are called. A single Gaussian is a **PoC simplification**. Multimodal haplotype structure is why the [evaluation](https://toki-bio.github.io/knowlege/projects/alu-graph-map-sv-proposal/) prefers mixture / haplotype-aware edges.

---

## The two formulas, and what they actually are

### Uniqueness weight \(W_u\)

For locus \(L\), tile canonical 31-mers; let \(C(k)\) be the genome-wide count of k-mer \(k\):

\[
W_u(L) = \frac{\lvert\{k \in K(L) : C(k) = 1\}\rvert}{\lvert K(L)\rvert}
\]

The spec then writes an IID misplacement probability for a read of length \(r\):

\[
P(\text{no unique }k\text{-mer}) = (1 - W_u)^{r-k+1}
\]

and claims \(W_u \ge 0.8\) is a **mathematical guarantee** that a 150 bp read cannot be misplaced (the formula at \(W_u=0.8\), \(k=31\) is astronomically small).

Treat that as a **bound under an assumption that Alu sequence violates**. K-mers inside an Alu body are highly shared and spatially correlated. A locus can have \(W_u = 0.8\) because unique k-mers sit in **diverged positions or immediate flanks**, while a read that lives entirely in the conserved core still multi-maps. \(W_u\) of the *annotated interval* is not uniqueness of *every overlapping read*.

What remains useful: \(W_u\) as a **ranking / filter** (old diverged AluJ vs young AluY). Expected “35–45% of Alus pass \(W_u \ge 0.8\)” is a prediction to measure with a k-mer database, not a theorem.

### “Topological entropy” \(H_\text{norm}\)

For each informative edge with baseline \((\mu_{ij},\sigma_{ij})\) and observed spacer length \(L_{ij}\):

\[
Z_{ij} = \frac{L_{ij} - \mu_{ij}}{\sigma_{ij}}, \qquad
H = \sum_{ij} Z_{ij}^2, \qquad
H_\text{norm} = H / N_\text{used}
\]

Edges whose endpoints fail a \(W_u\) threshold are excluded from \(H\). Design target: normal genomes \(H_\text{norm} \approx 1\); tumors \(\gg 1\); COLO829 tumor/normal fold \(\ge 2\) (hoped 4–8×).

This is a **sum of squared z-scores**, i.e. a chi-squared-like scatter statistic, not Shannon entropy and not a topological invariant of the graph. The name is branding.

If independent \(Z \sim N(0,1)\), then \(E[Z^2] = 1\) and \(H_\text{norm} \approx 1\). Edges are **not** independent (one SV hits many edges; mapping errors cluster; purity scales all Z). So “normal = 1” is a calibration target, not an identity. Tumor-vs-matched-normal fold is the honest first test; absolute \(H_\text{norm}\) across unrelated people is harder.

### How an SV is supposed to look

The spec’s local decoder:

```text
Z ≪ 0     spacer shortened     → deletion-like
Z ≫ 0     spacer lengthened    → duplication / insertion-like
|Z| ~ 0   within baseline
no spanning reads
where baseline expects them     → adjacency loss (e.g. translocation)
```

A **large** deletion that removes the whole spacer (or one of the Alus) will often look like **missing adjacency**, not a tidy \(Z = -5\). Small spacer indels can look like modest Z without being cancer SVs. Interpret Z as “this edge is weird,” then inspect reads — do not treat the cutoffs (−5 / +5) as biology.

Somatic **Alu insertions** into a spacer inflate \(L_\text{obs}\) (false +Z). The spec’s mitigation: look at spanning sequence for an Alu motif.

---

## Measurement stack (what actually observes an edge)

| Tier | What it sees | Honest limit |
| --- | --- | --- |
| Short reads (Illumina/BGI) | junctions / presence at **high-\(W_u\)** nodes | cannot span multi-kb spacers |
| Outward Alu-PCR | enriched short spacers | deletion-biased; Alu–Alu PCR chimeras |
| ONT low-pass 1–3× | molecule spanning node–spacer–node | the topology measurement the graph wants |
| ONT adaptive sampling (Read Until) | more on-target depth at ~WGS-run cost | keep/reject library must **not** require the healthy spacer (or novel junctions get rejected) |

**Sample type is a hard split**, not a footnote:

* **Solid tumor / long molecules** — ONT can span \(v_i \to e_{ij} \to v_j\). This is v1.
* **cfDNA (~160 bp)** — cannot span edges. Detect *discordant Alu pairs* with short reads / PCR, or drop the claim. Read Until does not rescue ultra-short fragments.

The COLO829 in-silico PoC (PacBio CCS) tests the **long-read spanning** path, not liquid biopsy.

---

## Pipeline (once, then per sample)

```text
once:
  Alu BED + genome k-mer DB  →  Wu per locus
  HPRC long-read haplotypes  →  baseline (μ, σ, n, orientation) per edge
                                (need ≥~10 haplotypes with spanning evidence)

per sample:
  spanning (or junction) evidence per edge
        ↓
  Z vs baseline; drop low-Wu endpoints from H
        ↓
  H_norm  +  per-edge outlier list  +  hotspot subgraph tests
```

Baseline construction on **short-read** HPRC BAMs will not populate 3–8 kb edges. Use haplotype-resolved long reads / pangenome paths. Sigma floor (e.g. 5 bp) is there so tiny technical jitter does not explode Z.

---

## COLO829 proof of concept (the first honest test)

**Data (as specified):** COLO829 (high-CIN melanoma) vs COLO829BL (matched lymphoblast), PacBio CCS, GRCh38, ENA PRJEB36890.

**Scripts, in order:**

1. `build_baseline.py` — HPRC haplotypes → `pangenome_baseline.tsv` (once).
2. `alu_wu_threshold.py` — 31-mer counts → \(W_u\) and the fraction of AluJ/S/Y above 0.8.
3. `alu_entropy_poc.py` — tumor + normal vs baseline → *H*_norm, Z histograms, overlap with published COLO829 SVs.

**PoC succeeds only if all of these hold**, not if a plot looks busy:

* \(H_\text{norm}(\text{tumor}) / H_\text{norm}(\text{normal}) \ge 2\)
* fraction of \(\lvert Z\rvert > 3\) clearly higher in tumor than normal
* DI/ID edges enriched among tumor outliers (NAHR-ish, not required for CIN in general)
* a large share of high-\(\lvert Z\rvert\) tumor edges sit near published COLO829 SV calls (the spec says ≥50% within ±50 kb — that window is wide; report a tighter overlap too)

If tumor and normal do not separate after dropping low-\(W_u\) nodes, the global score is not a biomarker yet. Local hotspot tests can still be worth a narrower product.

---

## Clinical outputs (intended)

1. **\(H_\text{norm}\)** — agnostic CIN / SV-load biomarker from low-pass data. Analogized to TMB, but for topology. Correlation with prognosis is **claimed**, not shown.
2. **Hotspot subgraphs** — posterior that a driver neighborhood’s edges sit outside the pangenome interval, with read support.
3. **Joint structural–epigenetic index** (optional ONT extension) — Alu hypomethylation plus \(H_\text{norm}\) from one run. Biologically motivated; extra claim until COLO829 / TCGA-style checks exist.

Subclonal SVs (low VAF) will not move bulk \(H_\text{norm}\) much. Position it as **macro-architecture screening**, not a resistance monitor.

Alu-sparse sequence (heterochromatin, some clinical loci) leaves long unanchored edges. Flag low-density zones; do not imply uniform SV sensitivity. L1/SVA nodes are a later densification, not a free upgrade.

HPRC ancestry is not a world census. Matched tumor/normal (COLO829) is robust to that; a universal \(H_\text{norm}\) threshold is not.

---

## Novelty vs US 10,242,154 — as claimed, not as legal advice

Closest prior art named in the spec: US 10,242,154, Alu/SVA/LINE markers for cancer, built around **inter-Alu PCR** and a **matched personal reference**, with inverted-Alu classifiers.

This spec’s own contrast:

| Piece | Older patent (as described) | This spec |
| --- | --- | --- |
| Markers | Alu (+ SVA/LINE) | Alu, divergence-weighted \(W_u\) |
| Reference | matched personal | population pangenome for the global score |
| PCR | core | one optional tier |
| Graph + pangenome priors | no | yes |
| Global CIN scatter score | no | central |
| ONT adaptive sampling / methylation joint index | no (era) | yes |

Four “pillars” listed for a possible filing: the Wu-weighted probabilistic graph; \(H_\text{norm}\) as SV-load biomarker; Read Until on spacers; joint methylation+structure index.

That is a **working comparison**, not a freedom-to-operate opinion. PCR overlap with the older patent is the part that needs a lawyer, not a slogan. The spec itself says: do not give a public talk or preprint on the assay until that process is decided; PoC evidence would be COLO829 in silico. None of that is a scientific result.

---

## Roadmap (only the scientific part)

```text
Phase 1  in silico     Wu on GRCh38; HPRC baseline; COLO829 fold + SV overlap
Phase 2  more lines    CIN rank order (e.g. HCC1187, K562) if Phase 1 separates
Phase 3  tissue        only after the statistic survives cell lines
```

Do not lead with 200-patient FFPE/cfDNA or FDA language. The evaluation page’s go/no-go list still applies: Read Until must keep *novel* junctions; entropy must not be chimera; sensitivity must be reported **by SV class and Alu density**.

---

## Bottom line

```text
landmarks = old, unique-enough Alus
signal    = spacer length / adjacency vs pangenome
score     = sum of squared Z (not entropy)
v1 data   = long molecules (COLO829 CCS), not cfDNA
Wu        = useful filter; the IID “guarantee” is not a guarantee
Gaussian  = PoC; mixtures later
patent    = process, not a result
```

Keep this page as the **spec**. Keep the [evaluation](https://toki-bio.github.io/knowlege/projects/alu-graph-map-sv-proposal/) as the **failure-mode list**. They should stay in disagreement where the spec is optimistic.
