#!/bin/sh
# Jake Yeung
# 4-rename_barcode_files.sh
# Rename barcode_counts.txt with the sample name 
# 2018-12-15

maindir="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/barcode_summaries_BM_build95"

for d in $(ls -d $maindir/PZ*); do
    dbase=$(basename $d)
    inf=$d/barcode_counts/bc_counts.txt
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    ln -s $inf $outdir/$dbase.bc_counts.txt
done
