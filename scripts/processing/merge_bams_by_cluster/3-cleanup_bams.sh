#!/bin/sh
# Jake Yeung
# 3-cleanup_bams.sh
#
# 2019-02-04

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged"

for c in `seq 1 8`; do
    bfile="$indir/H3K4me1_cluster_$c.bam"
    [[ ! -e $bfile ]] && echo "$bfile not found, exiting" && exit 1
    rm $bfile
done
