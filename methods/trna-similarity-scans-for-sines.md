---
title: Scanning SINEs for tRNA (and 5S / 7SL) similarity
type: explanation
---

> Type: explanation (deep research 2026-10-18, sourced) plus a section of measured results on rsi.

# Scanning SINE consensuses for tRNA, 5S and 7SL similarity

To decide whether a SINE consensus has a tRNA, 5S rRNA or 7SL RNA head, run three detectors that fail in different ways. The first is a **sequence-similarity test against a curated set of real Pol III RNA genes**: the SINEBase criterion is **at least 60% identity over at least 60 nt**, and the isotype is whichever reference RNA aligns best. The second is a **covariance-model scan** (tRNAscan-SE 2.0 or Rfam models in Infernal), read as **raw bit scores, not at default cutoffs**. The third is a **profile-HMM scan against Dfam's tRNA and SINE-head families**. Promoter boxes, folding and covariation only corroborate. None of them proves a head on its own, and none of them disproves one. The default cutoffs are where degraded heads get lost. They were set to keep SINEs *out* of tRNA gene sets: the tRNAscan-SE authors place tRNA-derived remnants in a **20–50 bit "twilight" zone**, and its high-confidence filter needs **70 bits on the isotype model**. The local rsi run shows this. Rfam and tRNAscan-SE missed the degenerate heads of r2–r4, and Dfam's tRNA-Ile model found them at E ≈ 1e-6 to 1e-7. The practical consequences follow. Run every detector with relaxed reporting and calibrate it against shuffled and genomic negative controls. Assign the isotype from whole-molecule models, not from the anticodon, and report the margin over the next-best isotype. When nothing is found, write "head origin not detected at threshold X", never "not tRNA-derived". No published benchmark measures any of these tools on degraded SINE heads, so every threshold below has to be checked against your own controls.

## The field's working test is 60% identity over 60 nt to a real RNA

The SINE-biology school (Kramerov, Vassetzky and colleagues) assigns heads by similarity, not by motifs. SINEBase admits a family only if it is an interspersed repeat with at least 100 copies and **"at least 60% identity with a tRNA species, 5S rRNA or 7SL RNA in at least 60-nt overlap"**. The exceptions are a few families whose Pol III transcription was confirmed experimentally. The same paper says A and B boxes "can serve only as an indication (but not a proof)" ([Vassetzky & Kramerov 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3531059)). Its search tool, SINESearch, runs FASTA against an RNA bank of human tRNAs, 5S and 7SL, and selects hits on **overlap length and identity rather than E-value**. The authors' reason is that E-values on short, degenerate similarities favour short near-perfect matches and miss real hits when the bank holds many similar sequences ([PMC3531059](https://pmc.ncbi.nlm.nih.gov/articles/PMC3531059)). Isotype is assigned the same way: Das-1's progenitor was named from its "significant similarity to alanine tRNA" ([Borodulina & Kramerov 2005](https://doi.org/10.1016/j.gene.2004.12.035)). Passerine heads were assigned to tRNA-Ile and tRNA-Asp by alignment, and cloverleaf drawings with shaded boxes were shown only as support ([Suh et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5729268)). By that count **84% of the 175 SINEBase families are tRNA-derived**, and **13% are composite** ([PMC3531059](https://pmc.ncbi.nlm.nih.gov/articles/PMC3531059)).

### The reference set decides what "best hit" means

The current tRNA reference is GtRNAdb **Release 22 (September 2024)** ([GtRNAdb](https://gtrnadb.ucsc.edu/)). Each per-genome tarball separates the high-confidence genes (`*-tRNAs.fa`, DNA, no CCA) from the pseudo and filtered predictions. The human hg38 tarball holds 619 predictions, of which 432 are in the confidence set and 187 are pseudo or filtered ([hg38 tarball](https://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz)). For bats, the only genome found was *Myotis lucifugus* (myoLuc2, 2010) ([Mluci2 summary](https://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Mluci2/Mluci2-summary.html)). This was not a full census; GtRNAdb may hold other bat genomes under names that were not searched. Two traps follow from the file layout. First, use only the high-confidence FASTA. Many of the pseudo and filtered loci are themselves SINE copies, so with them in the bank the "best tRNA hit" is just another SINE. Second, skip the "mature" FASTA, which is written in the U alphabet with an added 3′ CCA that genomic copies lack. The Leipzig tRNAdb that SINEBase was built on answered only with redirects when tested on 2026-10-18, and its status is listed as defunct ([Database Commons](https://ngdc.cncb.ac.cn:443/databasecommons/database/id/4011)). GtRNAdb is the practical replacement.

### 5S and 7SL heads need their own references and no B-box requirement

tRNA and 7SL heads carry a type-2 promoter (A box plus B box). 5S heads carry a type-1 promoter (A box, IE, C box) and **no B box**. A B-box-dependent scan therefore misses 5S SINEs, and bats have several of them. In megabats, MEG-RS and MEG-RL are 5S-derived, MEG-T2 is tRNA-derived and **MEG-TR is a tRNA/5S hybrid**. The 5S families keep their promoter better than their secondary structure ([Gogolevsky et al. 2009](https://doi.org/10.1016/j.ygeno.2009.02.001)). Repbase lists *Pteropus* Meg_PVa as SINE3/5S ([Repbase](https://girinst.org/2010/vol10/issue4/Meg_PVa.html)). Several other heads are composite. Tree-shrew Tu and bushbaby GarnAlu combine a tRNA part with a 7SL part, and hagfish SINE1-1_EBu combines a 7SL part with an older tRNA SINE ([Kojima 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7245038)). A fourth head class, SINEU, comes from U1/U2 snRNA ([Kojima 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4494074)). Because of these composites, a head test should score each RNA class along sliding windows of the first ~150–250 nt rather than keep the single best hit. On rsi, two of the elements turned out to be **tRNA–tRNA dimers**: r1_r3 has tRNA-Ile at 2–75 and again at 195–262, and r5h_r6 has it at 5–77 and 141–210 ([REFINEMENT.md §13](../../Tal/rsi/REFINEMENT.md)). A best-hit-only test reports only the first unit of each.

## Covariance models lose degraded heads at their default cutoffs

tRNAscan-SE 2.0 is the strongest structural tool, but its defaults are tuned *against* SINEs. It reports hits at **20 bits** by default and suggests 30 bits for large genomes. It calls a pseudogene when the total score is **below 55 bits** and either the primary score is below 10 or the secondary-structure score is below 5. Its eukaryotic high-confidence filter drops anything with an **isotype score below 70, a total below 50 or a structure score below 10** ([Chan et al. 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8450103/)). In practice, more than 80% of raw mammalian predictions are pseudogenes: cat has over 403,500 predictions against 549 high-confidence genes. The filter removed about 4,000 mouse repeats (mean 27.8 bits) and about 17,000 minke whale repeats (mean 48.0 bits). The authors also give the diagnostic: tRNA-derived SINEs keep the internal promoter but lose cloverleaf base pairing, and this shows up as a **drop in the secondary-structure score**. In the same paper, of 704 false positives found in shuffled "virtual" genomes at a 10-bit cutoff, **687 scored 10–20 bits** and the maximum was 29.5 bits. That gives the noise floor ([PMC8450103](https://pmc.ncbi.nlm.nih.gov/articles/PMC8450103/)). tRNAscan-SE can be relaxed: `-X` lowers the cutoff, `--max` removes the HMM filter (the web help recommends it for short, atypical sequences), `-H` splits the score into primary and structure parts, `-D` turns off pseudogene calling, and `-s` writes the isotype-specific model scan ([man page 2.0.9](https://hexmos.com/freedevtools/man-pages/user-commands/file-management/trnascan-se/)).

The Rfam models behave the same way. In Rfam 15.1, the gathering thresholds (GA) are **tRNA RF00005 29 bits (model length 71)**, **tRNA-Sec RF01852 47 bits**, **5S RF00001 38 bits** and **Metazoa_SRP/7SL RF00017 70 bits (model length 299)** ([RF00005 CM](https://rfam.org/family/RF00005/cm)). GA is curated so that no known false positive scores above it, which also means `--cut_ga` rejects most degraded heads. The 7SL model is a special case. It spans the full 299-nt RNA, and an Alu-like head covers only its Alu domain (human 7SL regions 1–63, 76–83 and 267–299) ([Kojima 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7245038)). Even a genuine 7SL head therefore falls below GA 70, so read it by E-value. For fragments, `cmsearch`/`cmscan --anytrunc` allows truncated hits that start and end anywhere, at the cost of less accurate E-values. `--max` disables truncation detection, so it has to be run as a separate pass ([cmsearch man page](https://dyn.manpages.debian.org/testing/infernal/cmsearch.1)). ARAGORN is a fast second opinion (99% sensitivity on 1,290 genes) with a single threshold knob (`-ps<num>`) and no isotype model ([Laslett & Canback 2004](https://pmc.ncbi.nlm.nih.gov/articles/PMC373265)).

### The rsi run shows the failure directly

The 2026-09-29 run compared three searches on the rsi consensuses: cmscan against all of Rfam (E ≤ 1e-3, relaxed to E ≤ 1), tRNAscan-SE `-E` at its defaults, and Dfam with `--cut_ga`. The controls came out as their names say: MEG-RS/RL as 5S, MEG-T2 as tRNA-Val, and MEG-TR as tRNA-Val (1–56) plus 5S (65–179). The rsi families split by how degraded their heads are ([REFINEMENT.md §13](../../Tal/rsi/REFINEMENT.md); [COMPOSITES.md §6](../../Surface/toki/SINEderella/docs/COMPOSITES.md)).

| Query | Dfam (tRNA-Ile model) | Rfam | tRNAscan-SE |
|---|---|---|---|
| r9 | 5–79, E 1e-23 | tRNA, E 3e-11 | Ile-AAT 65.4 bits, not pseudo |
| r7, r8 | 5–79, E 1e-18 to 1e-19 | tRNA | Lys / Met, pseudo |
| r5, r6 | 5–77, E 1e-10 to 1e-11 | none | Ser / Leu, pseudo |
| r1–r4 | 2–75, E 1e-6 to 1e-7 | none | r1 Ser (pseudo); **r2–r4 none** |
| r1_r3 linker (157–195) | none | none | none |

Three lessons follow from this table. First, the **detection order tracks degradation**. All three tools find r9. Rfam drops out at r5/r6, and tRNAscan-SE at its defaults drops out at r2–r4. Dfam, a sequence-profile HMM, still finds every head. One likely explanation, not yet tested: a covariance model rewards base pairs that old heads have lost, while a profile HMM scores only the conserved sequence, including the promoter boxes that survive. Second, the **anticodon calls are noise**. tRNAscan-SE read Ser, Lys, Met, Arg and Leu anticodons from heads that the Dfam model and the isotype models all place as Ile. The tRNAscan-SE authors call this kind of anticodon/isotype disagreement a "strong indicator" of decaying tRNA-derived repeats ([PMC8450103](https://pmc.ncbi.nlm.nih.gov/articles/PMC8450103/)). Third, the run had **controls on both sides**. Known positives (MEG families, r9) and a stretch of known-unrelated sequence (the r1_r3 linker, matched by nothing) bracket the result. That is the minimum calibration any head call needs.

## How the annotation tools label heads, read from code

Most SINE finders treat the head label as a tag, not a test. In AnnoSINE_v2 (cloned 2026-10-18), Step 4 runs `blastn -word_size 7 -evalue 1 -reward 2 -penalty -3 -gapopen 5 -gapextend 2` against a built-in bank of 1,729 sequences: 954 tRNA, 712 5S and 63 7SL. These counts equal the Rfam seed counts for RF00005 and RF00001. Each candidate gets the class of its **first** hit, or "Unknown" if there is none, and Unknown candidates are kept. Candidates with a best E-value **below 1e-15 are removed** as likely real RNA genes. The code comment says 1e-20, but the code uses 1e-15. The structure mode is SINE-Finder's regex: A box `[GA][CGA]TGG`, a 25–50 nt spacer, B box `GTTC[AG]A`, then a poly-A. These are 5–6-nt motif cores, so they match often by chance. The homology mode runs nhmmer with 118 Dfam SINE family HMMs, not with tRNA, 5S or 7SL models ([AnnoSINE_v2](https://github.com/liaoherui/AnnoSINE_v2); [Mobile DNA 2024](https://doi.org/10.1186/s13100-024-00331-y)). SINE_Scan uses an "enhanced SINE-Finder" that adds 7SL and 5S types to the original tRNA-only regex. It validates candidates by copy number and by a "similar SINE region flanked by unrelated sequence" test ([Mao & Wang 2017](https://academic.oup.com/bioinformatics/article/33/5/743/2687284)). RepeatMasker does not re-derive heads at all. Its SINE/tRNA, SINE/5S and SINE/7SL classes are the curators' family labels, and SINEBase warns that RepeatMasker "often leads to confusing misidentifications" because it ranks short high-identity fragments above long lower-identity matches ([PMC3531059](https://pmc.ncbi.nlm.nih.gov/articles/PMC3531059)).

Dfam 4.0 (May 2026, 30,621 curated families) is useful in two separate ways ([Dfam API](https://dfam.org/api/version)). Its `Pseudogene;RNA;tRNA` class contains per-isotype tRNA models: the tRNA-Ile model is the one that found the rsi heads. Its SINE classes (`SINE;tRNA_Promoter`, with 1,028 mammalian families including MIR and LFSINE_Vert, and `SINE;tRNA_and_7SL_RNA`, with 112 including B4 and ID_B1) answer a different question: "which known SINE does this head resemble". On rsi, r10's head matched MamSINE1 at 8–104 ([REFINEMENT.md §13](../../Tal/rsi/REFINEMENT.md)). Report the two as separate columns. A match to a known SINE is not evidence of tRNA origin by itself, because that SINE's own head classification is a curator's call. Dfam's documentation of how it assigns these subclasses was not found.

## Motifs, folding and covariation corroborate but cannot prove

The human type-2 promoter consensus is **box A `TRGYnnARBGG` at about +12 to +22 and box B `GGTTCRAnYCY` at about +54 to +64**, two 11-bp boxes 30–40 bp apart. The gap is wider when the tRNA has a long variable arm. 7SL, B1 and Alu share an atypical box A ([Tatosyan et al. 2020](https://www.mdpi.com/1422-0067/21/10/3706)). Another source gives the B box as `AGGTTCCAnnCC` ([Genome Res 14:1711](https://genome.cshlp.org/content/14/9/1711/T1.expansion.html)), so a report must say which definition it uses. Rough arithmetic for uniform base composition: the degenerate B box matches about once per 65,000 positions per strand, so a bare box hit genome-wide means little. A box counts only at the right position in a tRNA-length frame, at the right spacing from box A, and with a score calibrated on the target genome. In mammals CpG depletion lowers the count of `TCGA` cores, which makes an observed B box more informative there. The boxes have historical weight: Rhin-1 was discovered by PCR with a box-A primer ([Borodulina & Kramerov 2005](https://doi.org/10.1016/j.gene.2004.12.035)).

Single-sequence folding of a 70–90 nt head is weak evidence. MFE prediction gets only **40–70% of base pairs right** ([ViennaRNA RNAalifold tutorial](https://viennarna.readthedocs.io/en/latest/tutorial/RNAalifold.html)). The null must also preserve dinucleotides. Workman and Krogh showed that the apparent extra folding stability of mRNAs disappears against dinucleotide-matched shuffles, because stacking energies depend on neighbouring bases ([uShuffle, Jiang et al. 2008](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2375906/); [Clote et al. 2005](https://people.scs.carleton.ca/~kranakis/Papers/Zscores.pdf)). In old SINEs the cloverleaf is often gone: in the MIR/CORE literature, tRNA similarity "remains often limited to boxes A and B" ([PMC15861](https://pmc.ncbi.nlm.nih.gov/articles/PMC15861)). A failed fold therefore says little. Covariation with R-scape is the principled structural test, but it is poorly posed on SINE copies. R-scape tests whether paired-column covariation exceeds what phylogeny explains, and it reports statistical power, so a non-significant result does not rule out structure ([Rivas et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5554622)). Most copy-to-copy variation in a SINE family is neutral decay after insertion, with no compensatory changes. Only the source genes or subfamily masters were ever under structural selection. The better-posed test is therefore R-scape on **subfamily consensuses**, not on genome copies. No published study was found that has applied it to SINEs.

## Calibrate every score against two nulls

No published study was found that calibrates tRNA covariance-model scores on degraded heads, or on shuffled or genomic 90-nt windows. These numbers have to be measured. Two nulls are needed, because each fails differently. **Dinucleotide shuffles** of the query (`esl-shuffle -d`, or `ushuffle -k 2`) are the right null for folding energy and for CM structure scores ([esl-shuffle](https://www.mankier.com/1/esl-shuffle)). For a similarity score, however, shuffles are too easy, because they destroy every motif. **GC-matched random genomic windows** of the same length, taken from repeat-masked regions of the same genome, give the realistic null for "more tRNA-like than genome background". Turn each score into an empirical p-value, p = (1 + #null ≥ observed)/(1 + N). With N = 1,000 the smallest reportable p is about 0.001. Apply Benjamini–Hochberg across the heads tested. Some hairpins stay stable even after shuffling, a documented failure mode of z-scores for GC-rich or repetitive sequence ([bioRxiv 2023.11.30.569461](https://www.biorxiv.org/content/10.1101/2023.11.30.569461.full.pdf)). Always report GC content next to any MFE z-score. On the positive side, use known heads from the same genome or clade: r9 (clean tRNA-Ile), the MEG families (tRNA-Val and 5S) and Dfam MIR. For a negative control, use the r1_r3 linker or other parts that match nothing, plus Alu or B1 as a 7SL decoy for the tRNA test.

Two more checks belong in the calibration. First, true RNA genes and retro-pseudogenes pass every similarity test, so a head call needs SINE evidence alongside it: copy number, unique flanks, and **substitutions shared among copies but absent from the source tRNA**. SINEBase gives a single 51-nt human "B2" as a likely tRNA pseudogene, not a SINE ([PMC3531059](https://pmc.ncbi.nlm.nih.gov/articles/PMC3531059)). AnnoSINE's E < 1e-15 exclusion is a crude version of this check. Second, score heads against the mitochondrial tRNA models (`tRNAscan-SE -M mammal`, 22 isotype CMs) and flag any head that scores better there than against nuclear tRNAs as a possible NUMT. No mammalian SINE with a confirmed mt-tRNA head was found in the literature, so treat this as a precaution, not a known case.

## Thresholds to use, and what each means

| Detector | Report everything above | Tier "supported" | Tier "conservative" | Source |
|---|---|---|---|---|
| ssearch36 vs GtRNAdb/5S/7SL | any local hit | ≥60% id over ≥60 nt, beats genomic null p < 0.01 | same, plus margin over the runner-up isotype | [SINEBase](https://pmc.ncbi.nlm.nih.gov/articles/PMC3531059) |
| tRNAscan-SE 2.0 | `-X 0`, `--max`, `-D` | ≥20 bits and above the shuffle null | ≥30 bits | [Chan 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8450103/) |
| cmscan, Rfam tRNA/5S/7SL/Sec | `--anytrunc -E 10` | E ≤ 1e-3 | `--cut_ga` (29/38/70/47 bits) | [Rfam CM](https://rfam.org/family/RF00005/cm) |
| nhmmer, Dfam tRNA + SINE HMMs | `-E 10` | E ≤ 1e-5 and above the genomic null | `--cut_ga` | [rsi §13](../../Tal/rsi/REFINEMENT.md) |
| A/B box PWM with spacing | best pair in frame | p < 0.01 vs genomic windows | — | [Tatosyan 2020](https://www.mdpi.com/1422-0067/21/10/3706) |
| RNAfold MFE z vs dinucleotide shuffles | all | z ≤ −2 with GC reported | — | [Clote 2005](https://people.scs.carleton.ca/~kranakis/Papers/Zscores.pdf) |

The "supported" cutoffs are starting points that the controls must confirm, not published values. The E ≤ 1e-5 cutoff for the Dfam HMMs puts the weakest rsi heads (1e-6 to 1e-7) just inside it. Accept them only if no GC-matched genomic window scores as well. Call a head "tRNA-derived" when at least two independent detector families agree: sequence (ssearch36 or Dfam) and structure (a CM). Call it "tRNA-derived, sequence only" when only the sequence detectors pass. Call the isotype only when the isotype CM and the best GtRNAdb hit agree, or when one of them leads its runner-up by a margin that the shuffles rarely reach.

## A reproducible recipe for bat genomes

The script below assumes the therioserver environment already recorded in [COMPOSITES.md §6](../../Surface/toki/SINEderella/docs/COMPOSITES.md): conda `rnatools` with Infernal and tRNAscan-SE, and `~/refs/smallrna/Rfam.cm`. It adds ViennaRNA, HMMER, FASTA36 (ssearch36) and a local Dfam, so the run works offline instead of through the web API. Flags taken from the man pages cited above are verified. The `famdb.py` call and the R-scape conversion are written from working knowledge and should be checked against the installed versions; a toy run with the positive and negative controls must pass before any real run.

```bash
#!/usr/bin/env bash
set -euo pipefail
# --- inputs ----------------------------------------------------------------
W=/home/toki/heads_scan_$(date +%Y%m%d)   # therioserver: work under /home/toki, never /var/tmp or /tmp
mkdir -p "$W"; cd "$W" || exit 1
CONS=$1          # SINE consensuses (rsi bank) + controls (r9, MEG-RS, MEG-T2, MEG-TR, r1_r3 linker, AluY)
GENOME=$2        # repeat-masked Rhinolophus sinicus assembly (soft- or hard-masked)
REF=~/refs/smallrna; NSHUF=1000; THREADS=16
export TMPDIR=$W/tmp; mkdir -p "$TMPDIR"

# --- 0. references (once) --------------------------------------------------
# GtRNAdb Release 22: high-confidence DNA FASTA only (not *-filtered-*, not *-mature-*)
#   myoLuc2-tRNAs.tar.gz (bat), hg38-tRNAs.tar.gz (human); check counts vs summary page after tar tzf
cat $REF/myoLuc2-tRNAs.fa $REF/hg38-tRNAs.fa > $REF/trna_hc.fa
# 5S and 7SL genes (e.g. human RNA5S1, RN7SL1) + U1/U2 snRNA -> $REF/pol3_other.fa
cat $REF/trna_hc.fa $REF/pol3_other.fa > $REF/pol3_bank.fa
cmfetch -f $REF/Rfam.cm RF00005 RF01852 RF00001 RF00017 > $REF/pol3.cm  # one accession per line if -f wants a file
cmpress -F $REF/pol3.cm
# Dfam tRNA pseudogene models + mammalian SINE-head models (check famdb.py syntax for the installed version)
famdb.py -i $REF/Dfam.h5 families --format hmm --class "Pseudogene;RNA;tRNA" > $REF/dfam_trna.hmm
famdb.py -i $REF/Dfam.h5 families --format hmm --ancestors --descendants --class SINE Mammalia > $REF/dfam_sine.hmm

# --- 1. heads and sliding windows (composites: tRNA+tRNA, tRNA+5S) ---------
python3 - "$CONS" <<'PY' > windows.fa
import sys
seqs={};name=None
for l in open(sys.argv[1],encoding="utf-8"):
    l=l.strip()
    if l.startswith(">"): name=l[1:].split()[0]; seqs[name]=[]
    elif name: seqs[name].append(l.upper())
for n,s in seqs.items():
    s="".join(s)
    for st in range(0,max(1,len(s)-90+1),15):          # 90-nt windows, step 15, whole consensus
        print(f">{n}|{st+1}-{min(st+90,len(s))}\n{s[st:st+90]}")
PY

# --- 2. negative controls --------------------------------------------------
esl-shuffle -d -N $NSHUF windows.fa > null_shuf.fa                     # dinucleotide-preserving
# GC-matched (+/-2%) 90-nt windows from masked-free regions of the same genome (write with your sampler;
# exclude any window overlapping a lowercase/N base) -> null_gen.fa, NSHUF per GC bin

# --- 3. detectors (run identically on windows.fa, null_shuf.fa, null_gen.fa) ---------
for Q in windows null_shuf null_gen; do
  ssearch36 -n -3 -E 10 -m 8C -T $THREADS $Q.fa $REF/pol3_bank.fa > $Q.ssearch.m8   # -3: forward strand of query only
  cmscan --anytrunc -E 10 --cpu $THREADS --tblout $Q.cm.tbl $REF/pol3.cm $Q.fa > /dev/null
  cmscan --max      -E 10 --cpu $THREADS --tblout $Q.cmmax.tbl $REF/pol3.cm $Q.fa > /dev/null
  nhmmer -E 10 --cpu $THREADS --tblout $Q.dfam_trna.tbl $REF/dfam_trna.hmm $Q.fa > /dev/null
  nhmmer -E 10 --cpu $THREADS --tblout $Q.dfam_sine.tbl $REF/dfam_sine.hmm $Q.fa > /dev/null
done
# tRNAscan-SE on whole consensuses (it finds its own frame); relaxed + isotype scan + score breakdown
tRNAscan-SE -E -X 0 -H -D --max --thread $THREADS -o cons.trna.out -f cons.trna.ss -s cons.iso.txt "$CONS"
tRNAscan-SE -E -X 10 -H --max --thread $THREADS -o cons.trna_pseudo.out "$CONS"   # keeps pseudo labels
tRNAscan-SE -M mammal -X 0 -H --thread $THREADS -o cons.mt.out "$CONS"             # NUMT flag
tRNAscan-SE -E -X 0 -H -D --max --thread $THREADS -o shuf.trna.out null_shuf.fa    # null for bits

# --- 4. structure (corroboration only) -------------------------------------
RNAfold -p --noLP < windows.fa > windows.fold     # MFE per window; same on null_shuf.fa for z-scores
# Family level: align heads of SUBFAMILY CONSENSUSES (not genome copies), then
#   RNAalifold -p --aln heads.aln ; convert to Stockholm with SS_cons ; R-scape -s heads.sto

# --- 5. summarise: per window and per detector -> score, empirical p, best class/isotype, margin --
python3 summarise_heads.py --obs windows --nulls null_shuf null_gen \
        --trnascan cons.trna.out cons.iso.txt cons.mt.out --min-id 60 --min-ovl 60 > heads_report.tsv
```

The summariser (`summarise_heads.py`, still to be written) carries the logic of this report. For each window and detector, it takes the best hit per RNA class (tRNA isotype, 5S, 7SL, snRNA, Dfam SINE family). It computes the empirical p against both nulls, for that detector only, and stores the best-versus-second-best margin. It then joins overlapping windows into segments so that a consensus can come out as "tRNA-Ile 2–75 + tRNA-Ile 195–262" (r1_r3) or "tRNA-Val 1–56 + 5S 65–179" (MEG-TR). Each segment gets one of the tiers in the threshold table, together with the tRNAscan-SE primary and structure split and its IPD flag. From tRNAscan-SE, the summariser should read the isotype CM name and score, not the anticodon. A head that scores better on `-M mammal` than on `-E` should be flagged as a possible NUMT. A consensus with no segment above its null must be written out as "head origin not detected at threshold X", followed by the detector list. Before any real run, the toy run must reproduce the known rsi table above: all r-families tRNA-Ile, MEG-RS 5S, MEG-TR tRNA plus 5S, and nothing on the r1_r3 linker. Only then is the pipeline trustworthy on a new bat genome.

## Conclusion

The rsi episode shows that "Rfam and tRNAscan-SE found nothing" says more about thresholds and model type than about the sequence. The default cutoffs were built to purge SINEs from tRNA gene catalogues, and structure-rewarding models lose old heads before sequence-profile models do. A head scan should therefore lead with sensitive sequence detectors (pairwise identity against GtRNAdb, Dfam HMMs) and use covariance models and folding as a second, independent check. It should name isotypes from whole-molecule models with an explicit margin, and treat anticodon calls as untrustworthy for anything older than r9.

The largest open gap is not a missing tool but a missing calibration. Nobody has published detector sensitivity as a function of head divergence. Bat genomes offer a natural benchmark: a graded series from clean (r9) to barely detectable (r2–r4), plus 5S and composite heads (MEG-RS, MEG-TR). Measuring every detector against these heads and against the same two nulls would turn the tiers above from rules of thumb into measured operating points, and that benchmark could itself be a methods result worth reporting.


## Applied: rsi (Rhinolophus sinicus), 2026-10-18 (measured)

42 queries (21 consensuses + the stretches of 7 composites). Files and full table: Tal
`rsi_main_20261018/trna/` (README, `trna_summary.tsv`, script `trna_rsi.sh`).

- ssearch36 against the 419 human tRNA genes (GtRNAdb hg19) found a **tRNA-Ile head in every r-unit**: 73-93 % over
  66-75 nt, bit score 39-82 against at most 29 for 20 dinucleotide shuffles of the same query; next amino acid 8-37 bits
  lower. All pass SINEBase's 60 % / 60 nt rule.
- **Rfam missed r2, r3, r4, r6; tRNAscan-SE missed r3, r4** - the pairwise search against tRNA genes, with the
  shuffle null, found all of them. Anticodon calls of the decayed heads disagree with the isotype model (IPD).
- Composites carry one tRNA head per unit (P1: two; C11: three; P26 / P34 / P48: two; P18: one, its middle and its r8-like
  part have none). MEG-RL / RS = 5S, MEG-T2 = tRNA-Val, MEG-TR = tRNA-Val + 5S, as published.
