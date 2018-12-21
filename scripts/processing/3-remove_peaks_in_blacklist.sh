#!/bin/sh
# Jake Yeung
# 3-remove_peaks_in_blacklist.sh
# Remove peaks from blacklist region 
# 2018-12-18

blacklist="/hpc/hub_oudenaarden/jyeung/data/scChiC/blacklist/mm10.blacklist.bed.gz"
# inbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_macs2_output/BM_H3K4me1_merged.0.3.1000.cutoff/BM_H3K4me1_merged.0.3.1000.cutoff_peaks.broadPeak"
# outbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_macs2_output/BM_H3K4me1_merged.0.3.1000.cutoff/BM_H3K4me1_merged.0.3.1000.cutoff_peaks.blacklistfilt.broadPeak"
# do 100kb merge to compare with MetaCell

dist=25000
# inbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_macs2_output/BM_H3K4me1_merged.0.3.1000.cutoff/BM_H3K4me1_merged.0.3.1000.cutoff_peaks.blacklistfilt.broadPeak"
inbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_macs2_output/BM_H3K4me1_merged.0.3.1000.cutoff/BM_H3K4me1_merged.0.3.1000.cutoff_peaks.merge_${dist}.broadPeak"
outbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_macs2_output/BM_H3K4me1_merged.0.3.1000.cutoff/BM_H3K4me1_merged.0.3.1000.cutoff_peaks.merge_${dist}.blacklistfilt.broadPeak"

[[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1
[[ ! -e $blacklist ]] && echo "$blacklist not found, exiting" && exit 1
[[ -e $outbed ]] && echo "$outbed not found, not overwriting so exit" && exit 1

bedtools intersect -v -a $inbed -b $blacklist > $outbed
