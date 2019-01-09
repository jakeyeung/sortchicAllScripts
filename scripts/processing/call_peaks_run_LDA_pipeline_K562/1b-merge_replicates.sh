#!/bin/sh
# Jake Yeung
# 1b-merge_replicates.sh
# Top replicates together. Can do very stringent filtering for calling peaks 
# When actually analyzing data we can use a less stringent filteing 
# 2019-01-08

dmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_K562/all_bam_links"

[[ ! -d $dmain ]] && echo "$dmain not found, exiting" && exit 1

# only two marks are done apparently
marks="H3K27me3 H3K4me3"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged_K562"

[[ ! -d $outdir ]] && mkdir $outdir

for mark in $marks; do
    echo $mark
    inbams=$(echo $dmain/PZ-K562-$mark-*.bam)
    outbam=$outdir/K562_"$mark"_merged.bam
    [[ -e $outbam ]] && echo "$outbam exists, skipping" && continue
    samtools merge $outbam $inbams
    samtools index $outbam
done




