#!/bin/sh
# Jake Yeung
# 0-rename_bigwigs.sh
# Rename neutrophil bigwigs 
# 2019-03-21

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Neutrophils/Gong_GenesAndDev_2017/bigwig_mm10"
outdir=$indir/renamed

for b in `ls -d $indir/*.bw`; do
    bname=$(basename $b)
    bname=${bname%.*}
    # keep mark and celltype only
    mark=$(echo $bname | cut -d"_" -f2)
    celltype=$(echo $bname | cut -d"_" -f3)
    bnew=${mark}_${celltype}
    bout=$outdir/$bnew.bw
    ln -s $b $bout
done
