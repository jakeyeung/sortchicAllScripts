#!/bin/sh
# Jake Yeung
# 2a-remove_header_from_good_cells.sh
# Remove header  
# 2019-05-09

pcutoffcell=0.05
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/quality_control/good_cells.pcounts.${pcutoffcell}.pfrac.${pcutoffcell}.txt"
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/quality_control/good_cells.pcounts.${pcutoffcell}.pfrac.${pcutoffcell}.noheader.txt"

[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

awk 'NR>1' $inf > $outf
