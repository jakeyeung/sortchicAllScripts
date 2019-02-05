#!/bin/sh
# Jake Yeung
# 2019-02-04_copy_merged_bams.sh
# Visualize merged bams initially 
# 2019-02-04

indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged"
outdir="/Users/yeung/data/scchic/bam_cluster_merged"

rsync -avrL --copy-links $indir $outdir
