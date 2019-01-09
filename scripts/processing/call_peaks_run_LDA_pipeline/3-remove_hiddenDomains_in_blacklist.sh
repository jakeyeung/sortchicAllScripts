#!/bin/sh
# Jake Yeung
# 3-remove_hiddenDomains_in_blacklist.sh
# Remove peaks from blacklist region 
# 2019-01-08


marks="H3K4me1 H3K27me3 H3K9me3 H3K4me3"
mindist="1000"
cell="BM"

for jchip in $marks; do
    echo $jchip
    inbeddir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_bam_hiddenDomains_output"
    inbedname="${cell}_${jchip}_merged.${mindist}.cutoff/${cell}_${jchip}_merged.${mindist}.cutoff_analysis.bed"
    inbed=$inbeddir/$inbedname
    blacklist="/hpc/hub_oudenaarden/jyeung/data/scChiC/blacklist/mm10.blacklist.bed.gz"
    outbedname="${cell}_${jchip}_merged.${mindist}.cutoff/${cell}_${jchip}_merged.${mindist}.cutoff_analysis.blacklistfilt.bed"
    outbed=$inbeddir/$outbedname

    [[ -e $outbed ]] && echo "$outbed found, skipping" && continue
    [[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1
    [[ ! -e $blacklist ]] && echo "$blacklist not found, exiting" && exit 1
    # [[ -e $outbed ]] && echo "$outbed found, not overwriting so exit" && exit 1
    bedtools intersect -v -a $inbed -b $blacklist > $outbed
done
