#!/bin/sh
# Jake Yeung
# 1d-rename_tagged_dir.sh
# rename tagged dir 
# 2019-05-10

indir="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95-B6"
n=0
maxjobs=16
for d in `ls -d $indir/B6*`; do
    bname=$(basename $d)
    indir=$d/tagged
    outdir=$d/tagged_old2
    mv $indir $outdir
done
wait

