# Pol III A/B box retention in SINE consensi (polIII_boxes_20261025)

Server work directory (actual run, 2026-10-25): KIT:/data/V/toki/polIII_boxes_20261025/ (therioserver and DRAGEN were unreachable that day;
run_polIII_boxes.sh is the untested therioserver variant, run_polIII_boxes_kit.sh is what ran)
Local copy of code and outputs: C:\work\Raks_COI\research_notes\Pol III A and B boxes in SINEs\data\

## Question
How well do SINE consensus sequences (SINEBase) keep the Pol III type-2 promoter boxes (A box, tRNA
positions 8-18; B box, positions 52-62), compared with real human tRNA genes and with a
dinucleotide-shuffled null?

## Files
| file | what |
|---|---|
| clean_fasta.py | first header token as name, uppercase, U->T, gaps removed |
| box_pwm.py | PWM building from tRNA genes, scoring of consensi and shuffles, all tables |
| run_polIII_boxes_kit.sh | the pipeline that ran on KIT (inputs -> 3 Rfam CMs -> esl-shuffle -> cmscan -> box_pwm.py; introns from GtRNAdb hg19-tRNAs-detailed.out) |
| run_polIII_boxes.sh | therioserver variant (full Rfam.cm, tRNAscan-SE) - written first, NOT run |
| inputs/ | SINEBase.nr95.uniq.fa (local SINEBase copy, 229 consensi), GtRNAdb hg19 tarball (hg19-tRNAs.fa 419 genes, hg19-tRNAs-detailed.out), rfam/RF00005+RF00001+RF00017, md5.txt |
| out/per_consensus.tsv | one row per consensus (scores, positions, 11-mers, Galli mismatches, flags, dimer, head columns) |
| out/summary_by_head.tsv | counts and % by Rfam head type; last row = tRNA genes |
| out/A_without_B.tsv, out/B_without_A.tsv, out/dimers.tsv, out/controls.tsv | subsets of per_consensus.tsv |
| out/trna_boxes.tsv | per tRNA gene: chosen A/B sites, leave-one-out scores, spacing, intron |
| out/trna_summary.txt | PWM consensus/IC, tRNA score percentiles, spacing distribution, null percentiles, flag definitions |
| out/pwm_A.tsv, out/pwm_B.tsv | the two log-odds matrices (pseudocount 0.5, uniform background) |
| out/cm.tbl | cmscan --tblout vs tRNA, 5S_rRNA, Metazoa_SRP models (E <= 1) |
| out/sines.shuf.fa | esl-shuffle -d -N 20 --seed 42 |

## Commands
```
# KIT, through the one ssh tunnel (plink -P 2223 toki@127.0.0.1):
pscp clean_fasta.py box_pwm.py run_polIII_boxes_kit.sh README.md SINEBase.nr95.uniq.fa -> /data/V/toki/polIII_boxes_20261025/
cd /data/V/toki/polIII_boxes_20261025 && bash run_polIII_boxes_kit.sh > run.log 2>&1
```

## Method in one paragraph
A and B box sites are seeded in each tRNA gene by the loose regexes T[AG]G[CT]..A[AG].GG (window: start
4-30) and GGTTC[AG]A..C[CT] (window: start within 34..18 nt of the gene end); genes without a regex hit
get the 11-mer with the fewest mismatches to the Galli consensus (TGGCNNAGTGG / GGTTCGANNCC) in the
window. Log-odds PWMs (pseudocount 0.5, background 0.25) are built, each site is re-located once with the
PWM inside its window, and the PWMs are rebuilt. tRNA-gene scores are leave-one-out. Consensi are scanned
on the sense strand: best A anywhere, best B anywhere, best B 20-120 nt downstream of the best A
(spacing), the same for the best A within the first 40 nt (head columns), and a second A/B pair downstream
of the first (dimers). Null: 20 dinucleotide-preserving shuffles per consensus (esl-shuffle -d, seed 42, Easel 0.48),
best A/B score per shuffle. "Present" = above the 95th percentile of the 20 shuffles of that consensus
(*_null) or at or above the 5th percentile of the tRNA-gene site scores (*_tRNA).

Column head_SINEBase in per_consensus.tsv is left empty for the SINEBase website family annotations (collected separately).
