#!/bin/sh
# Jake Yeung
# 2019-02-17_copy_bams_merged.sh
# Copy bams to explore
# 2019-02-17

inmain="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged"
outmain="/Users/yeung/data/scchic/bams_merged/."

scp -r $inmain $outmain/
