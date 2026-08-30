---
title: Tuatara tRNA-derived minisatellite project
type: notes
---

> Type: notes

# Tuatara tRNA-derived minisatellite — project reference

This is a **project record**, not an established result. Observations, hypotheses, goals, and planned analyses are kept separate on purpose so a working hypothesis is not mistaken for a conclusion months later.

---

## 1. Central observation

A class of **minisatellite / tandem-repeat sequences** has been identified in the tuatara genome using tandem-repeat detection.

The sequences appear **potentially** related to tRNA-derived sequence. Two non-mutually-exclusive possibilities remain open:

1. the minisatellite originated from **tRNA pseudogenes / tRNA-derived genomic fragments**, which later became tandemly amplified;
2. the minisatellite originated from an **ancient tRNA-derived SINE or another SINE-like element**, with tandem structure representing later amplification / evolution.

**Origin is not established.** The project exists to test the evolutionary relationship, not to assume tRNA origin from superficial similarity.

---

## 2. Biological question

> What is the evolutionary origin of this tuatara minisatellite family, and what is its relationship to tRNA genes, tRNA pseudogenes, and ancient SINEs?

More specifically:

* Are the repeat monomers genuinely derived from tRNA sequence?
* Does the entire monomer correspond to a recognizable tRNA-derived structure, or only part of it?
* Are different repeat copies descendants of one ancestral tRNA-derived element?
* Could the repeats instead be remnants of an ancient SINE that was itself tRNA-derived?
* Is the family unique to tuatara, or are there homologs among reptiles / birds?
* What does its distribution imply about age and evolutionary history?

---

## 3. Existing tuatara data

### Genome

```text
/data/W/toki/Genomes/Reptiles/Tuatara/GCA_003113815.1_ASM311381v1_genomic.fna
```

Assembly:

```text
GCA_003113815.1_ASM311381v1
```

### TRF output

```text
trf GCA_003113815.1_ASM311381v1_genomic.fna 2 5 7 80 10 50 2000 -h
```

Primary `.dat` file:

```text
GCA_003113815.1_ASM311381v1_genomic.fna.2.5.7.80.10.50.2000.dat
```

This `.dat` file is currently the primary source for genomic repeat loci.

---

## 4. Preliminary consensus

A preliminary consensus exists, but it is explicitly **not a final family consensus** and has **not been demultiplexed**.

```text
>Consensus
AGGGTCCAGGGTTCAAGTCCCTGTGCAGGCGGCNGGTCTCCCTCHHTCTTGCATGCTCCCTAGGCCAACTTTTGTGATGGTCTAGCGGAGGTTTAGAGAGAAACAAACAAACAAACAAGAAAAGACACCCTGAGGAANACTGCAAGCAGAGCCTGGGTAGCTCAGCCGGTAGAGCATCAGACTCTTAATCTGAGGGTCCAGGGTTCAAGTCCCTGTTCAGGCATCAGCTCTCCTGGGAGTTGGACTCTCCAACNCACATCCCAGGHCTNGCTTCTGTGGGGGGAGGTGAGTTGGCTGTGTTTTTGCCTTCACCGGTGCCCACCACCTCNTCCCCCCHCCCCCCN
```

Characteristics to remember:

* substantial sequence complexity;
* ambiguous bases (`N`, `H`, etc.);
* apparent low-complexity / poly-A-like region in the middle;
* another GC-rich / low-complexity region toward the end;
* length looks compatible with a tRNA-derived element **plus additional sequence**, not necessarily a clean canonical tRNA;
* that interpretation must be tested, not inferred from appearance.

Treat this consensus as a **starting clue**, not as the biological family’s sequence.

---

## 5. Planned overall strategy

The intended analysis is modeled on the earlier **CLsat** work, not on generic repeat annotation alone.

```text
TRF calls
   ↓
filter appropriate monomer / period length
   ↓
identify repeat arrays / loci
   ↓
demultiplex related repeat types / subfamilies
   ↓
derive proper consensus sequences
   ↓
estimate genomic abundance
   ↓
extend genomic loci
   ↓
analyze extended loci and flanking sequence
   ↓
compare repeat copies / consensuses with tRNAs
   ↓
compare across species
   ↓
compare with other repetitive elements / SINEs
   ↓
infer evolutionary origin and history
```

Principle:

> Characterize the repeat family on its own terms first.  
> Only then force a tRNA-derived or SINE-derived interpretation.

This is not simply “CLsat with another repeat.” The interesting biological possibility is an evolutionary connection of the form:

**tRNA → repetitive / mobile DNA → minisatellite**

---

## 6. Stage A — define the repeat family correctly

First practical task: inspect the TRF output and determine the distribution of **period / monomer lengths**.

Questions:

* Is there a clear dominant period?
* Is there more than one period peak?
* Are apparent different periods variants of the same underlying repeat?
* Does TRF identify one family or several related structures?
* Are some TRF calls fragments of larger arrays?
* What thresholds separate genuine members from incidental tandem repeats?

### Filtering

Filter the `.dat` file by the monomer length / period characteristic of this family, in the same conceptual manner used for CLsat.

Other TRF parameters (copy number, score, array length, etc.) should be considered only after the period distribution is understood.

Parsers must follow the **actual TRF `.dat` format**, not assumed generic column positions.

---

## 7. Locus definition

After selecting relevant TRF hits, define genomic **repeat loci / arrays**.

Important distinctions:

* one TRF hit ≠ one biological locus;
* adjacent / overlapping TRF calls may be parts of one tandem array;
* a single locus may contain many imperfect monomers;
* interrupted arrays may produce multiple TRF calls.

Document a reproducible definition of:

> What constitutes one minisatellite locus?

Those loci become the units for abundance, distribution, extraction, consensus construction, extension, and between-locus comparisons.

---

## 8. Demultiplexing — unresolved and important

The preliminary consensus was produced **without demultiplexing**.

Different related variants may currently be collapsed into one sequence.

After locus extraction, examine whether there are:

* one homogeneous family;
* several closely related subfamilies;
* different evolutionary generations;
* structurally distinct variants;
* a conserved ancestral component plus rapidly evolving additions / deletions.

What looks like one family may resolve into:

```text
                 ┌── subfamily A
ancestral repeat ├── subfamily B
                 └── subfamily C
```

or into structurally different classes. That can change the whole origin story.

---

## 9. Consensus reconstruction

Once family / subfamilies are defined, build proper consensuses.

Goals:

* which positions are conserved;
* which regions are variable;
* whether variation concentrates in particular parts;
* whether the apparent poly-A / low-complexity segment is genuine;
* whether conserved blocks correspond to portions of tRNA;
* whether subfamilies preserve different parts of a putative ancestral sequence.

Desired products:

1. a family consensus;
2. potentially several subfamily consensuses;
3. alignments of genomic copies showing variation.

---

## 10. Abundance analysis

Once loci are defined, estimate abundance in the tuatara genome:

* number of loci;
* number of arrays;
* number of repeat units / copies;
* total genomic bp occupied;
* fraction of the genome occupied;
* array-length distribution;
* monomer-number distribution;
* scaffold / genomic distribution.

Distinguish:

**locus abundance** vs **repeat-copy abundance**.

A few huge arrays and thousands of short arrays are different evolutionary histories even if total copy number looks similar.

---

## 11. Extension of loci

As in the CLsat analysis, extend identified loci into their genomic surroundings.

Flanks may contain the best origin evidence:

* remnants of flanking SINE / LINE sequence;
* target-site duplications;
* poly-A tails;
* insertion boundaries;
* neighboring repeats;
* degradation patterns;
* evidence of an ancestral mobile-element insertion;
* whether tandem amplification occurred inside a pre-existing repeat context.

Extended sequences are not just extra alignment material. Their **genomic architecture** needs examination.

---

## 12. tRNA comparison

Primary comparative hypothesis: tRNA origin.

Compare against **actual tRNA genes**, not only generic tRNA consensus sequences.

Initial reference scope:

### Tuatara

* genomic tRNA genes;
* identifiable tRNA pseudogenes, if available.

### Other reptiles

* squamates;
* turtles;
* crocodilians.

### Birds

Useful as an additional archosaur comparison outside reptiles sensu stricto.

Exact species set can be chosen later by genome quality and availability.

---

## 13. What tRNA comparison must establish

A single BLAST hit is not enough.

Need:

### Sequence-level consistency

Do independent minisatellite copies match the same tRNA-derived regions?

### Structural consistency

Does the monomer correspond to recognizable tRNA architecture — for example acceptor stem, D arm/loop, anticodon region, variable region, T arm, or other conserved components?

### Evolutionary consistency

Do relationships among repeat copies and tRNAs make biological sense?

Desired model:

```text
tRNA
  ↓
ancestral tRNA-derived element
  ↓
amplification
  ↓
minisatellite family
```

rather than only:

```text
unrelated repetitive sequence
       ↕
generic tRNA similarity
```

---

## 14. Ancient SINE hypothesis

Investigate independently — complementary, not automatic.

Questions:

* recognizable SINE architecture?
* evidence of a tRNA-derived SINE lineage?
* partner LINE / retrotransposon-associated fragments?
* poly-A tails?
* target-site duplications?
* insertion boundaries in extended loci?
* related but non-tandem copies elsewhere in the genome?
* some “minisatellite units” actually remnants of an ancestral mobile element?

Keep open:

> **tRNA → SINE → tandem amplification**

versus only:

> **tRNA → tandem minisatellite**

These are different biological hypotheses and need different genomic evidence.

---

## 15. Cross-species analysis

Comparative scope:

```text
Tuatara
   │
   ├── Squamates
   ├── Turtles
   ├── Crocodilians
   └── Birds
```

Ask more than “does an identical repeat exist?”

Investigate:

1. exact / similar repeat sequences;
2. tRNA-derived sequence similarity;
3. homologous repeat families;
4. related SINEs;
5. conserved fragments of an ancestral element;
6. genomic distribution;
7. phylogenetic presence / absence pattern.

A highly divergent homolog elsewhere may be more informative than absence of a near-identical repeat.

---

## 16. Other repetitive elements

Also compare against broader repeat datasets:

* SINEs;
* LINE-associated sequences;
* known tRNA-derived retrotransposons;
* other minisatellites;
* elements with similar architecture.

Candidate outcomes to keep distinct:

```text
tRNA-derived minisatellite
SINE-derived minisatellite
tRNA-derived SINE → minisatellite
another evolutionary scenario
```

---

## 17. Most important conceptual constraint

Do **not** begin by assuming the preliminary consensus is the ancestral sequence.

Instead:

```text
observed genomic copies
        ↓
family structure
        ↓
subfamilies
        ↓
consensus(s)
        ↓
ancestral relationships
        ↓
comparison to tRNA / SINE candidates
```

Tandem repeats evolve by substitutions, indels, unequal crossing-over, slippage, homogenization, expansion / contraction, and sequence turnover. A present-day consensus can differ substantially from the sequence that founded the family.

---

## 18. Current status

### Already available

* tuatara genome;
* TRF `.dat` output;
* preliminary repeat consensus;
* working hypotheses about tRNA / pseudogene / SINE origin;
* CLsat-inspired analytical strategy.

### Not yet established

* exact monomer period;
* optimal TRF filtering criteria;
* number of loci;
* abundance;
* number of subfamilies;
* proper consensus;
* locus-extension statistics;
* relationship to tuatara tRNAs;
* relationship to tRNAs in other taxa;
* SINE relationship;
* evolutionary age / distribution;
* final origin model.

---

## 19. Recommended order of work

```text
A. Inspect TRF output
B. Determine period distribution
C. Define / filter candidate repeat hits
D. Define genomic loci
E. Characterize locus and array statistics
F. Extract individual repeat arrays
G. Demultiplex into families / subfamilies
H. Build proper consensuses
I. Extend representative / all loci
J. Analyze flanking genomic context
K. Identify tuatara tRNAs
L. Compare minisatellite ↔ tRNAs
M. Search for SINE / mobile-element relationships
N. Comparative analysis in other taxa
O. Build evolutionary model
```

Do not jump directly to phylogeny or declare “tRNA-derived.” The first job is to establish exactly what repeat family is under study.

---

## 20. Immediate next step

Start from the **actual TRF `.dat` file**, not from a generic pipeline.

Inspect representative `Sequence:` blocks and TRF repeat records and establish:

1. exact TRF field structure;
2. period distribution;
3. candidate period range;
4. whether multiple repeat classes are present;
5. how much of the preliminary consensus is actually represented across independent loci.

Only after that write extraction / filtering code. That keeps the project from becoming a pile of arbitrary thresholds and scripts.
