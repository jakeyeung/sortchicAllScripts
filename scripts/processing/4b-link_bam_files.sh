#!/bin/sh
# Jake Yeung
# 4b-link_bam_files.sh
# Create symlinks to bam files
# 2018-12-15

inmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam"

for d in $(ls -d $inmain/PZ*); do
    echo $d
    dbase=$(basename $d)
    bamf=$d/$dbase.filtered.sorted.bam
    bamfi=$d/$dbase.filtered.sorted.bam.bai
    [[ ! -e $bamf ]] && echo "$bamf not found, exiting" && exit 1
    [[ ! -e $bamfi ]] && echo "$bamfi not found, exiting" && exit 1
    ln -s $bamf $outmain/$dbase.bam
    ln -s $bamfi $outmain/$dbase.bam.bai
done
