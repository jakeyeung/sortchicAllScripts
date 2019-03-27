#!/bin/sh
# Jake Yeung
# 1b-merge_replicates.sh
# Top replicates together. Can do very stringent filtering for calling peaks 
# When actually analyzing data we can use a less stringent filteing 
# 2018-12-18

# dmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods"
dmain="/hpc/hub_oudenaarden/jyeung/data/histone-mods-Ensembl95"
# mark="H3K4me1"
# mark="H3K27me3"
marks="H3K4me1 H3K27me3 H3K9me3 H3K4me3"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged_build95"

[[ ! -d $dmain ]] && echo "$dmain not found, exiting" && exit 1
[[ ! -d $outmain ]] && mkdir $outmain

for mark in $marks; do
    echo $mark
    inbams=$(echo $dmain/PZ*$mark*/*sorted.bam)
    outbam=${outmain}/BM_"$mark"_merged.bam
    [[ -e $outbam ]] && echo "$outbam exists, skipping" && continue
    samtools merge $outbam $inbams
done




