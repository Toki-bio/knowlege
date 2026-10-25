#!/bin/bash
# run_polIII_boxes_kit.sh - Pol III A/B box retention in SINEBase consensi vs human tRNA genes vs dinucleotide-shuffled null.
# KIT version (therioserver and DRAGEN were unreachable on 2026-10-25). Work dir /data/V/toki/polIII_boxes_20261025/.
# Inputs : inputs/SINEBase.nr95.uniq.fa (local SINEBase copy, 229 consensi, transferred by pscp),
#          inputs/hg19-tRNAs.fa + hg19-tRNAs-detailed.out (GtRNAdb hg19 tarball, wget on KIT; 419 genes, introns kept),
#          inputs/rfam/RF00005 RF00001 RF00017 covariance models (tRNA, 5S_rRNA, Metazoa_SRP) fetched from rfam.org.
# Tools  : /data/V/toki/miniforge3/bin/python3 (3.10), /usr/local/bin/esl-shuffle, cmscan/cmpress from env stold.
# Usage  : cd /data/V/toki/polIII_boxes_20261025 && bash run_polIII_boxes_kit.sh > run.log 2>&1
set -euo pipefail
W=/data/V/toki/polIII_boxes_20261025
cd "$W" || exit 1
export TMPDIR=$W/tmp TMP=$W/tmp TEMP=$W/tmp
mkdir -p "$TMPDIR" inputs inputs/rfam out
PY=/data/V/toki/miniforge3/bin/python3
CMSCAN=/data/V/toki/miniforge3/envs/stold/bin/cmscan
CMPRESS=/data/V/toki/miniforge3/envs/stold/bin/cmpress
ESL=/usr/local/bin/esl-shuffle
echo "=== start $(date)  host $(hostname)  dir $W"
$PY --version; $ESL -h | head -2; $CMSCAN -h | head -2

# 1. inputs
[ -s inputs/hg19-tRNAs.fa ] || { cd inputs && wget -q -T 60 -O hg19-tRNAs.tar.gz "http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/hg19-tRNAs.tar.gz" && tar xzf hg19-tRNAs.tar.gz; cd "$W"; }
$PY clean_fasta.py inputs/SINEBase.nr95.uniq.fa inputs/sines.fa
echo "sines: $(grep -c '>' inputs/sines.fa)   tRNAs: $(grep -c '>' inputs/hg19-tRNAs.fa)"
md5sum inputs/SINEBase.nr95.uniq.fa inputs/hg19-tRNAs.fa inputs/hg19-tRNAs-detailed.out > inputs/md5.txt
for acc in RF00005 RF00001 RF00017; do
  [ -s inputs/rfam/$acc.cm ] || curl -sS -L -m 180 -o inputs/rfam/$acc.cm "https://rfam.org/family/$acc/cm"
done
cat inputs/rfam/RF00005.cm inputs/rfam/RF00001.cm inputs/rfam/RF00017.cm > inputs/rfam/rfam3.cm
echo "Rfam models: $(grep -c '^INFERNAL' inputs/rfam/rfam3.cm) ; names: $(grep '^NAME' inputs/rfam/rfam3.cm | awk '{print $2}' | sort -u | tr '\n' ' ')"
$CMPRESS -F inputs/rfam/rfam3.cm > /dev/null

# 2. null: 20 dinucleotide-preserving shuffles per consensus, fixed seed
$ESL -d -N 20 --seed 42 inputs/sines.fa > out/sines.shuf.fa
echo "shuffles: $(grep -c '>' out/sines.shuf.fa)"

# 3. head labels: cmscan vs the three models, E <= 1
$CMSCAN --cpu 8 -E 1 --tblout out/cm.tbl -o /dev/null inputs/rfam/rfam3.cm inputs/sines.fa
echo "cmscan hits: $(grep -vc '^#' out/cm.tbl)"

# 4. PWMs, scoring, tables
$PY box_pwm.py --trna inputs/hg19-tRNAs.fa --trna_detailed inputs/hg19-tRNAs-detailed.out \
  --sines inputs/sines.fa --shuf out/sines.shuf.fa --cm out/cm.tbl --out out | tee out/box_pwm.stdout
echo "=== done $(date)"
ls -la out/
