#!/bin/sh
# Jake Yeung
# copy_bin_files_from_avo.sh
# COpy bin files from AvO 
# 2019-03-22

marks="H3K4me1 H3K4me3 H3K9me3"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_tables_raw"

for mark in $marks; do
    echo $mark
    inf="/hpc/hub_oudenaarden/avo/scChiC/metacell/BM-$mark-build95.txt"
    cp $inf $outdir/.
done
