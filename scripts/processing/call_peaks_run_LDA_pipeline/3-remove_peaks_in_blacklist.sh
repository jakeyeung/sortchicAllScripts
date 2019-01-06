#!/bin/sh
# Jake Yeung
# 3-remove_peaks_in_blacklist.sh
# Remove peaks from blacklist region 
# 2018-12-18


marks="H3K4me1 H3K27me3 H3K9me3 H3K4me3"
pvalcutoff="0.5"
mindist="1000"

for jchip in $marks; do
    echo $jchip
    inbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_macs2_output/BM_${jchip}_merged.${pvalcutoff}.${mindist}.cutoff/BM_${jchip}_merged.${pvalcutoff}.${mindist}.cutoff_peaks.broadPeak"
    blacklist="/hpc/hub_oudenaarden/jyeung/data/scChiC/blacklist/mm10.blacklist.bed.gz"
    outbed="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_macs2_output/BM_${jchip}_merged.${pvalcutoff}.${mindist}.cutoff/BM_${jchip}_merged.${pvalcutoff}.${mindist}.cutoff_peaks.blacklistfilt.broadPeak"

    [[ -e $outbed ]] && echo "$outbed found, skipping" && continue
    [[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1
    [[ ! -e $blacklist ]] && echo "$blacklist not found, exiting" && exit 1
    # [[ -e $outbed ]] && echo "$outbed found, not overwriting so exit" && exit 1
    bedtools intersect -v -a $inbed -b $blacklist > $outbed
done
