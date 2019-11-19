#!/bin/sh
# Jake Yeung
# 2019-06-25_copy_H3K4me3_bams.sh
# Copy bams clustered already  
# 2019-06-25

indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bigwigs_build95_B6_stringent"
outdir="/Users/yeung/data/scchic/for_episys"

scp -r $indir $outdir
