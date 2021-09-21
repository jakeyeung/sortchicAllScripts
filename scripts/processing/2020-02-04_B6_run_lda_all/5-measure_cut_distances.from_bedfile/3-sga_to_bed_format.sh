#!/bin/sh
# Jake Yeung
# 3-sga_to_bed_format.sh
#  
# 2020-08-06

inf="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.sga.chromorenamed.gz"
outf="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTss.chromorenamed.bed.gz"

zcat $inf | awk 'BEGIN{OFS="\t"}{print ($1, $3, ($3 + 1), $4, $6)}' | gzip > $outf
# zcat $inf | cut -f1,3,4,6
# zcat $inf | cut -f1,3,4,5 | gzip > $outf
