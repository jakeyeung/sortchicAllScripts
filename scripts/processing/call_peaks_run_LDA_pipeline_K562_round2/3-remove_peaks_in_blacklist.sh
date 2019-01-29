#!/bin/sh
# Jake Yeung
# 3-remove_peaks_in_blacklist.sh
# Remove peaks from blacklist region 
# 2018-12-18


marks="H3K27me3 H3K4me3"
pvalcutoff="0.5"
mindist="1000"

cell="K562"
macsmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/macs2_output_K562"
for jchip in $marks; do
    echo $jchip
    beddir="$macsmain/${cell}_${jchip}_merged.${pvalcutoff}.${mindist}.cutoff"
    bedname="${cell}_${jchip}_merged.${pvalcutoff}.${mindist}.cutoff_peaks.broadPeak"
    bednameout="${cell}_${jchip}_merged.${pvalcutoff}.${mindist}.cutoff_peaks.blacklistfilt.broadPeak"
    inbed="$beddir/$bedname"
    blacklist="/hpc/hub_oudenaarden/jyeung/data/scChiC/blacklist/mm10.blacklist.bed.gz"
    outbed="$beddir/$bednameout"

    [[ -e $outbed ]] && echo "$outbed found, skipping" && continue
    [[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1
    [[ ! -e $blacklist ]] && echo "$blacklist not found, exiting" && exit 1
    # [[ -e $outbed ]] && echo "$outbed found, not overwriting so exit" && exit 1
    bedtools intersect -v -a $inbed -b $blacklist > $outbed
done
