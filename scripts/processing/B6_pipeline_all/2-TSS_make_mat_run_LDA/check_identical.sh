#!/bin/sh
# Jake Yeung
# check_identical.sh
#  
# 2019-05-09

indir="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6-bamlist"
# check are identical
marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

for mark in $marks; do
    fname="JY_${mark}_bamlist_goodcellsfilt.out"
    f1=$indir/$fname
    f2=$indir/windowbug/$fname
    diff $f1 $f2
done
