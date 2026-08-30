---
title: Crocodilian SINEs and MITEs
type: explanation
---

> Type: explanation

# SINEs in crocodile genomes — and why MITEs matter beside them

## Short answer

Yes, there is literature touching SINEs in crocodilians — but **few papers focus specifically on crocodilian SINE families**. Most crocodile genome work emphasizes overall TE landscapes dominated by ancient **CR1 LINEs**. SINEs are present but comparatively under-discussed: often ancient, low-activity, and annotation-poor.

MITEs (miniature inverted-repeat transposable elements) are a useful contrast: they can burst to high copy number in crocodilians (e.g. *Chompy*), even when SINE expansions look weak.

---

## What the crocodilian TE literature usually says

### Broad TE landscape

Crocodilian genomes are TE-rich (often cited around ~37% TE-derived sequence in vertebrate TE reviews / genome papers). The landscape is typically:

* dominated by **ancient elements**, especially **CR1 LINEs**;
* relatively **low recent TE activity**;
* fairly **conserved across crocodilian lineages**.

Working interpretation (useful, not dogma):

```text
TE “fossilization” environment
  → SINEs exist
  → but often ancient / degraded / lineage-specific
  → not primate-style Alu bursts
```

### CR1-centered amniote history

Papers on ancient CR1 lineages across crocodiles, turtles, and birds emphasize long-term stability of TE landscapes. Why that matters for SINEs:

> SINEs are non-autonomous and depend on LINE machinery.  
> If CR1 activity is low/stable, recent SINE mobilization is also likely suppressed.

### Functional TE fragments

Work on TE-derived microRNAs in saltwater crocodile (*Crocodylus porosus*) shows TE sequence can retain regulatory roles even when mobilization is quiet. That does **not** prove abundant active SINEs; it only shows TE debris can still matter functionally.

---

## Honest state of the field for SINEs specifically

### Known / expected

* SINEs are widespread across vertebrates (often tRNA- or 7SL-derived, Pol III, non-autonomous).
* Crocodilian genomes contain many ancient TE relics.
* TE activity is generally lower than in many mammals or some squamates.

### Sparse / missing

* No famous crocodile-specific SINE radiation analogous to primate *Alu* or rodent B1/B2.
* Few papers that:
  * classify crocodilian SINE families carefully,
  * quantify copy numbers by family,
  * reconstruct their evolutionary history in detail.

Likely biological + practical reasons:

* long generation times / relatively stable genomes,
* ancient degraded copies that are hard to annotate,
* study designs that report “% TE” bins rather than fine SINE taxonomy.

### Practical data sources if you want real numbers

Useful assemblies commonly include:

* *Alligator mississippiensis*
* *Crocodylus porosus*
* *Gavialis gangeticus*

Starting points:

* RepeatMasker / Dfam / Repbase annotations filtered for `SINE`,
* then your own [QuickSearch-style consensus search](https://toki-bio.github.io/knowlege/sines/sine-search-pipeline-pseudocode/) if curated libraries are incomplete.

Expect mostly ancient MIR-like / tRNA-derived material and lower recent expansion than LINE content suggests at first glance.

---

## Contrast case: MITEs in crocodilians

The *Chompy* story in alligator is informative even if it is not a SINE paper:

* a MITE-like family can reach tens of thousands of copies (~46k reported for *Chompy*);
* non-autonomous DNA elements **can** expand in crocodilians;
* by analogy, SINE expansions are biologically possible — they are simply less documented / apparently less dramatic in current annotations.

Useful contrast:

| Feature | Crocodilian SINEs (typical picture) | Crocodilian MITEs (e.g. Chompy) |
| --- | --- | --- |
| Recent expansion | weak / underreported | bursts possible |
| Dependency | LINE reverse transcriptase | DNA transposase |
| Annotation | often buried in ancient TE noise | structurally detectable (TIR + TSD) if you look |
| Signal type | evolutionary background | recent activity candidate |

---

## How MITEs multiply

MITEs encode no proteins. They hitchhike on autonomous DNA transposons.

### Mechanism

1. A full-length DNA transposon (e.g. hAT, Tc1/Mariner) encodes transposase.
2. Transposase recognizes **terminal inverted repeats (TIRs)**.
3. MITEs have similar TIRs but lack the internal coding region.
4. Transposase binds the MITE → excises → inserts elsewhere.

Classic DNA transposition is cut-and-paste, yet MITEs can reach huge copy numbers because host biology can preserve the donor copy (below).

### Hallmarks in data

* TIRs (~10–30 bp, often conserved),
* target-site duplications (TSDs; family-specific, e.g. TA for many Tc1/Mariner),
* short length (~100–600 bp),
* high copy number in bursty families.

### Why we know they are DNA elements

Multiple independent lines of evidence:

1. **TIRs** — retroelements like SINEs do not have this architecture.
2. **TSDs** matching DNA-transposon family signatures.
3. **Terminal similarity** to autonomous DNA transposons, with deleted internals.
4. **Mobilization experiments** — supply active transposase, MITEs move.
5. **No RNA-intermediate machinery** — no RT / typical SINE poly-A / Pol III package.

| Feature | MITEs | SINEs |
| --- | --- | --- |
| Mobility | DNA cut/paste (+ host repair) | RNA copy/paste |
| Requires | transposase | reverse transcriptase (LINE partner) |
| Ends | TIRs | no TIRs |
| TSDs | yes | yes (different insertion biology) |
| Expansion | bursts possible | depends on LINE activity |

---

## How can cut-and-paste raise copy number?

Naive puzzle:

```text
excise from A → insert at B
copy number unchanged?
```

Reality: donor sites often persist. Several routes have experimental support in DNA-transposon systems (Ac/Ds, P elements, yeast/plant/fly work). There is **not** one universal MITE-only law; several host processes can contribute.

### 1. Transposition during DNA replication

If excision happens after a region has been replicated, and insertion occurs into unreplicated DNA (or the complementary sister outcome retains a copy), net copy number can increase.

### 2. Gap repair after excision

Excision leaves a double-strand break. Homologous repair from a sister chromatid can restore the donor sequence while the excised copy inserts elsewhere.

### 3. Sister-chromatid / paired-template interactions

Related to the above: paired chromatids enable restoration + new insertion patterns and non-Mendelian outcomes.

What is **not** claimed: a single proven “MITEs always duplicate via X” rule. Amplification is the observed outcome; replication timing and repair are the best-supported explanations for preserving donor loci.

Practical inference when you see 10k–50k MITE copies:

* efficient mobilization,
* frequent donor restoration / replication-coupled gain,
* likely activity in dividing cells / S-phase contexts,
* then ask whether divergence profiles look like a burst or slow drip.

---

## Why so little MITE literature in squamates and tuatara?

The gap is real. It is mostly **annotation + study-design**, not proof of biological absence.

### What papers usually report

Squamate TE papers (including work arguing squamates challenge simple repeat-evolution paradigms) emphasize:

* highly variable TE content across species,
* LINE-dominated or lineage-specific mosaics (CR1, BovB, L2, etc.),
* high diversity with often modest copy number per family (“many families, each small”).

The tuatara genome paper (Gemmell et al.) reports a huge (~5 Gb) genome with a mixed reptile/mammal-like TE landscape, including notable SINE content — but again, **not a dedicated MITE catalog**.

### Why MITEs disappear from the narrative

1. **Annotation bias** — RepeatMasker/Repbase are stronger on retroelements; short fast-evolving MITEs get labeled `DNA/unknown`, fragmented, or missed.
2. **Study design** — papers ask “what % is TE?” and bin `DNA transposon`, so MITEs vanish inside a broad category.
3. **Different TE ecology** — crocs: fewer families, visible bursts; squamates: many families, each small → MITE signal diluted.
4. **DNA turnover** — mutation/deletion can erase TIRs and make old MITEs unrecognizable.
5. **Method mismatch** — homology annotation ≠ TIR/TSD structural discovery (MITE-Hunter / detectMITE / custom TIR scans).
6. **Tuatara weirdness** — giant genome, slow evolution, nested ancient debris; signal may be present but deeply diverged. See also the [tuatara minisatellite project notes](https://toki-bio.github.io/knowlege/projects/tuatara-trna-minisatellite/).

### Important counterpoint

DNA transposons are clearly present/active in squamate history (e.g. widely discussed SPIN-like elements). Where autonomous DNA transposons exist, MITE derivatives are expected. So:

```text
biological presence ≠ annotation visibility
```

This is a research gap more than a demonstrated absence.

---

## Implications for analysis

If mining reptile genomes:

* **SINEs in crocs** → expect ancient background; do not assume mammal-like expansions; verify with libraries + your own searches.
* **MITEs** → structurally detectable recent-activity candidates if you look for TIR + TSD + short high-copy clusters.
* **Squamates / tuatara** → absence of MITE papers is not absence of MITEs; structural discovery may reveal lineage-specific families.

A useful mental split:

```text
SINEs  → LINE-dependent retroecology
MITEs  → DNA-transposase + host repair ecology
```

---

## Key takeaways

* Dedicated crocodilian SINE literature is thin; TE papers mostly highlight ancient CR1-dominated landscapes.
* SINEs likely exist in crocs as old / low-activity / under-annotated elements rather than famous recent radiations.
* MITEs show that non-autonomous bursts *can* happen in crocodilians (*Chompy*).
* MITE copy-number growth is compatible with cut-and-paste via replication timing and gap repair that preserve donor sites.
* Sparse MITE reporting in squamates/tuatara is largely an annotation and research-focus gap, not proof they are missing.
* For your work: treat croc SINEs as a careful extraction problem, and treat MITEs as a structural-discovery opportunity beside LINE/SINE macro-stats.
