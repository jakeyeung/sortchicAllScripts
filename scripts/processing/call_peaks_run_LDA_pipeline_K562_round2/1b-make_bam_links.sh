#!/bin/sh
# Jake Yeung
# 1b-make_bam_links.sh
# Make bam links 
# 2019-01-15

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_K562_round2"
outdir="$inmain/all_bam_links"

[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
[[ ! -d $outdir ]] && mkdir $outdir

for f in `ls -d $inmain/PZ*/*sorted.bam*`; do
    echo $f
    ln -s $f $outdir/.
done

