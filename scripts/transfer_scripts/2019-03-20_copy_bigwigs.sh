#!/bin/sh
# Jake Yeung
# 2019-03-20_copy_bigwigs.sh
# Copy bigwigs  
# 2019-03-20

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bigwigs_2019-03-20"
outdir="/Users/yeung/data/scchic/from_cluster/bigwigs_2019-03-20"

scp -r t2:$indir $outdir
