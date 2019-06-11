#!/bin/sh
# Jake Yeung
# 2-create_grep_file_from_good_cells.sh
# Create grep patterns from good cells 
# 2019-05-08

inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/quality_control/good_cells.pcounts.0.1.pfrac.0.1.noheader.txt"
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/quality_control/good_cells.pcounts.0.1.pfrac.0.1.noheader.forgrep.txt"

# marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
# for mark in $marks; do
#     echo $mark
#     grep $mark $inf
# done
suffix=".sorted.bam$"
awk -v suffix=$suffix 'NF{print $0 suffix}' $inf > $outf
