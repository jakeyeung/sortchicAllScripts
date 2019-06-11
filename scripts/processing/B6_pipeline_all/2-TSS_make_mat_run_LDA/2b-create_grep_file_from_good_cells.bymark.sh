#!/bin/sh
# Jake Yeung
# 2-create_grep_file_from_good_cells.sh
# Create grep patterns from good cells 
# 2019-05-08

pcutoffcell=0.05
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/quality_control/good_cells.pcounts.${pcutoffcell}.pfrac.${pcutoffcell}.noheader.txt"

[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

suffix=".sorted.bam$"
marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
for mark in $marks; do
    outfmark="/hpc/hub_oudenaarden/jyeung/data/scChiC/quality_control/good_cells.pcounts.${pcutoffcell}.pfrac.${pcutoffcell}.noheader.forgrep.${mark}.txt"
    grep $mark $inf | awk -v suffix=$suffix '{print $0 suffix}' > $outfmark
done
# awk -v suffix=$suffix 'NF{print $0 suffix}' $inf > $outf
