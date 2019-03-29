#!/bin/sh
# Jake Yeung
# 2019-03-20_copy_clustered_bams.sh
# Copy clustered bams and explore 
# 2019-03-20

indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/sorted_bams_2019-03-20"
outdir="/Users/yeung/data/scchic/from_cluster"

scp -r $indir $outdir

