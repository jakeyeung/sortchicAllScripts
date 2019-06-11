#!/bin/sh
# Jake Yeung
# 3-filter_bamlist_by_good_cells.sh
# Filter by good cells 
# 2019-05-08

pcutoff=0.05
inmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6-bamlist"
# goodcells="/hpc/hub_oudenaarden/jyeung/data/scChiC/quality_control/good_cells.pcounts.0.1.pfrac.0.1.noheader.forgrep.txt"
goodcellmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/quality_control"

# [[ ! -e $goodcells ]] && echo "$goodcells not found, exiting" && exit 1
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

for mark in $marks; do
    fname="good_cells.pcounts.${pcutoff}.pfrac.${pcutoff}.noheader.forgrep.${mark}.txt"
    goodcells=$goodcellmain/$fname
    [[ ! -e $goodcells ]] && echo "$goodcells not found, exiting" && exit 1
    bamlist=$inmain/JY_${mark}_bamlist_all.out
    [[ ! -e $bamlist ]] && echo "$bamlist not found, exiting" && exit 1
    bamlistfilt=$inmain/JY_${mark}_bamlist_goodcellsfilt.out
    grep -f $goodcells $bamlist > $bamlistfilt&
done
wait
