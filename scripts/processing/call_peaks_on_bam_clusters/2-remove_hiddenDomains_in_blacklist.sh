#!/bin/sh
# Jake Yeung
# 2-remove_hiddenDomains_in_blacklist.sh
# Remove peaks from blacklist region 
# 2019-04-15


marks="H3K4me1 H3K27me3 H3K9me3 H3K4me3"
mindist="1000"
cell="BM"

inbeddir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/merged_cluster_bam_hiddenDomains_output_build95"
[[ ! -d $inbeddir ]] && echo "$inbeddir not found, exiting" && exit 1
blacklist="/hpc/hub_oudenaarden/jyeung/data/scChiC/blacklist/mm10.blacklist.bed.gz"
[[ ! -e $blacklist ]] && echo "$blacklist not found, exiting" && exit 1

for jchip in $marks; do
    echo $jchip
    for indir in $(ls -d $inbeddir/${jchip}_cluster_*.${mindist}.cutoff); do
        indirbase=$(basename $indir)
        inbedname=${indirbase}_analysis.bed
        inbed=$indir/$inbedname

        outbedname="${indirbase}"_analysis.blacklistfilt.bed
        outbed=$indir/$outbedname

        [[ -e $outbed ]] && echo "$outbed found, skipping" && continue
        [[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1
        [[ ! -e $blacklist ]] && echo "$blacklist not found, exiting" && exit 1
        # [[ -e $outbed ]] && echo "$outbed found, not overwriting so exit" && exit 1
        bedtools intersect -v -a $inbed -b $blacklist > $outbed

    done
done
