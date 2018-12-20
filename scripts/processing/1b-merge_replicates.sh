#!/bin/sh
# Jake Yeung
# 1b-merge_replicates.sh
# Top replicates together. Can do very stringent filtering for calling peaks 
# When actually analyzing data we can use a less stringent filteing 
# 2018-12-18

dmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods"
mark="H3K4me1"

[[ ! -d $dmain ]] && echo "$dmain not found, exiting" && exit 1

inbams=$(echo $dmain/PZ*$mark*/*sorted.bam)
# echo $inbams

outbam=/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged/BM_"$mark"_merged.bam
samtools merge $outbam $inbams



