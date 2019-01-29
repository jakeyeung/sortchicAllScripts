#!/bin/sh
# Jake Yeung
# 2d-merge_peaks_by_distance.sh
# Merge broad peaks called by Macs2 by distance for input into MetaCell 
# 2018-12-21

# set 50kb ends up being 100kb?
# dist=50000  # merge 100kb like AvO
dist=10000  # merge 100kb like AvO
# dist=25000  # merge 100kb like AvO

inbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_macs2_output/BM_H3K4me1_merged.0.3.1000.cutoff/BM_H3K4me1_merged.0.3.1000.cutoff_peaks.blacklistfilt.broadPeak"
outbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_macs2_output/BM_H3K4me1_merged.0.3.1000.cutoff/BM_H3K4me1_merged.0.3.1000.cutoff_peaks.merge_${dist}.broadPeak"

bedtools merge -d $dist -i $inbed -c 4,5,6,7,8,9 -o distinct,distinct,distinct,mean,mean,mean | awk -F $'\t' 'BEGIN {OFS = FS} $4=("Peak" FNR)'  > $outbed

echo "outpath: $outbed"
