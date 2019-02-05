#!/bin/sh
# Jake Yeung
# 2019-02-04_copy_hiddenDomains_bed.sh
#  
# 2019-02-04

jmark="H3K4me1"

inf="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_hiddenDomains_output/BM_${jmark}_merged.1000.cutoff/BM_${jmark}_merged.1000.cutoff_analysis.blacklistfilt.bed"

outf="/Users/yeung/data/scchic/beds/."

scp $inf $outf
