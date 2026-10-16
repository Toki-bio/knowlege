---
title: "tandem-array-unit-check — digest"
type: reference
topic: Program digests
generated: true
---

> Type: reference · Claude Code skill `tandem-array-unit-check` · [source on GitHub](https://github.com/Toki-bio/knowlege/tree/main/.claude/skills/tandem-array-unit-check)

**When to use:** How to decide cheaply whether a group of regularly spaced repeat copies (SINE hits, satellite monomers) is a real tandem array or chance clustering, without long pairwise alignments. Covers short flank / inside-piece tests, k-mer containment (Mash, ModDotPlot, AniAnn's), the limits of TRF, ULTRA, SRF, TRASH, minimap2 and megablast for kilobase units, and SINEderella's measured data (rsi, rle, Timema). Use when designing or reviewing any array / satellite / segmental-duplication check, or when a step aligns multi-kb stretches.

---

{% raw %}
# Is this run of copies a tandem array?

> Type: explanation (facts with sources) plus marked notes. Written 2026-10-13.

## The question

A run is >= 5 copies of one element at regular spacing. It is either:
- **a tandem array:** the stretch from one copy to the next is a copy of the same unit, so units are near-identical
  (rsi MEG-RS 2,155 bp unit, 99.9 %), or
- **chance:** ordinary dispersed copies that happen to sit near each other (unit identity 25–55 %, often mixed strands).

In hit-rich families chance runs vastly outnumber arrays. Timema t3-1: 1,996 runs checked, 12 arrays. rsi r-families:
1,850 of 1,862 five-to-nine-unit runs are copies.

## Owner's rule

No step compares long sequences (user, 2026-10-13). SINEderella's kind-B check aligned whole units (median 5 kb, up to
24 kb on Timema) with Smith-Waterman. That cost 85 of the 95 min of its satellite screen, so it is the example of what not
to do.

## Methods, cheapest first

| Method | What it compares | Cost | Caveat |
|---|---|---|---|
| **Flank check** | ~150 bp after (and before) each copy vs the next copy's and the one after (lag 2 for dimeric arrays) | 150 x 150 bp per pair | a copy at the very end of the unit has the next unit's copy as its flank; take both sides |
| **Inside piece** | ~150 bp from the middle of the inter-copy stretch, searched in the next stretch | 150 bp x unit length | tests that the stretch itself repeats |
| **k-mer containment** | share c of unit B's k-mers (both strands) found in unit A; identity ~ 100 x (1 + ln c / k) (Mash-style) | linear in length, no alignment | **fails near 70–85 %**: see below. Use only as a fast positive |
| minimap2 pair (`-c`, asm20) | seed-chain-align two units | ms per pair | seeds vanish at high divergence like k-mers; ModDotPlot authors found minimap2 slow on highly repetitive input |
| megablast | 28-mer words | per-call overhead | less sensitive near 85 %, no advantage |
| full SW (ssearch36) | whole units | grows with length x length | what to avoid; `-z -1` is 10x faster than `-z 11` but still quadratic |

## Measured facts (SINEderella `docs/SATELLITES.md` 5h–5i)

- **k-mer pre-screen rejected:** on 8,005 real rsi unit pairs (therioserver `~/tmp/kbspeed/calib.tsv`), 31 of 625 pairs scoring
  >= 70 had < 1 % of positions in shared 16-mers, and 7 shared no 12-mer. At 70–85 % identity the mismatches are 3–7 bp apart,
  so no 12–20-mer survives. The uniform-mutation expectation (85 % gives ~9 % shared 15-mers) does not hold on real old arrays.
- **TRF cannot do kb units:** TRF 4.10 refuses a maximum period above 2,000 bp; toy arrays of 2,500 and 7,000 bp units give no
  record.
- **Dimeric arrays exist:** rsi MEG-RS alternates 2,155 and 1,460 bp units sharing ~800 bp. Compare lag 1 **and** lag 2.
- **Run length is not the criterion, sequence is:** rle MEG-RS arrays are missed by the family share, and a run-length rule
  would take 1,850 rsi clusters of copies.
- **ssearch36 path limit:** it cannot open files whose path is longer than ~120 characters. Run inside a temp dir with bare names.

## Tools found (2026-10-13 search) and fit

- **AniAnn's** (2025): self-similarity matrix of ANI between genomic windows, k-mer containment sketched with modimizers, banded
  near the diagonal so O(n). Detects arrays from 32 bp units to rDNA (~45 kb). A whole-sequence scanner, not a candidate checker,
  but the same idea. https://pmc.ncbi.nlm.nih.gov/articles/PMC12873891/
- **ModDotPlot** (2024): interactive dotplot from k-mer containment; A. thaliana 128 Mb in < 5 min on a laptop.
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11071298/
- **StainedGlass**: minimap2 identity heatmaps; > 1 h for one 3 Mb centromere on a cluster.
  https://www.biorxiv.org/content/10.1101/2021.08.19.457003.full.pdf
- **Mash** distance D = -(1/k) ln(2j/(1+j)). https://mash.readthedocs.io/en/latest/distances.html
- **SRF** (Li 2023): de novo satellite units from a k-mer de Bruijn graph (a perfect array is one cycle). Reconstructs units;
  it does not judge given candidates. https://genome.cshlp.org/content/33/11/1994
- **TRASH** (2023): de novo tandem repeats + higher-order structure; aligns to a generated consensus, struggles with divergent
  and composite arrays (per AniAnn's). https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10199239/
- **ULTRA** (2024): HMM with one state per period and a maximum period; decayed repeats well, slow at kb periods.
  https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11580682/

## Recommended design (notes, not yet tested)

1. Verdict from **short pieces**: flank check (lag 1 + 2, both sides). Median per-unit score >= 85 gives ARRAY.
2. Inside-piece search as the second opinion when the flank check is borderline.
3. k-mer containment only as a **fast positive** (clearly ARRAY), never to reject, unless calibrated at the cut on real pairs.
4. Implement in bash/gawk + `samtools faidx -r` + short `ssearch36 -m 8`; **gawk, not mawk** (mawk is quadratic on long lines).
5. Validate verdict-by-verdict against the existing ssearch36 verdicts: rsi `validate_rsi.tsv` (47 ARRAY / 2,692 COPIES),
   rle MEG-RS/MEG-TR, Timema `satellites/loci.bed`. Toy-test first.

gawk containment sketch (lag 1; lag 2 = keep two previous units):

```awk
# input: run_id <TAB> unit_index <TAB> sequence, sorted
function rc(s,  i, r) { r = ""; for (i = length(s); i > 0; i--) r = r comp[substr(s, i, 1)]; return r }
BEGIN { FS = "\t"; k = 15; comp["A"]="T"; comp["T"]="A"; comp["C"]="G"; comp["G"]="C"; comp["N"]="N" }
{ s = toupper($3); n = length(s) - k + 1
  if ($1 == prun && n > 0) { hit = 0; r = rc(s)
    for (i = 1; i <= n; i++) if ((substr(s, i, k) in prev) || (substr(r, i, k) in prev)) hit++
    c = hit / n; print $1, $2 - 1, $2, c, (c > 0 ? 100 * (1 + log(c) / k) : 0) }
  delete prev; for (i = 1; i <= n; i++) prev[substr(s, i, k)]; prun = $1 }
```

## Worked example (Timema, consensus t3-1, run of 2026-10-13)

- ARRAY 99 %: `CM115228.1:114591405-114596530`, 7 hits, all `+`, spacing ~730 bp (SINE ~290 bp + ~440 bp).
- COPIES 8 %: `CM115225.1:7300916-7324215`, 6 hits, strands `+ - - - + +`, gaps 3.4–5.8 kb.
- Kind-B cost there: ~85 min for 14 consensuses; units median 5.1 kb, max 23.8 kb. The run is on KIT
  `/data/V/toki/tool_update_20261013/timema_new/run_20261013_090058`.
{% endraw %}
