#!/bin/sh
# Jake Yeung
# 1c-merge_replicates.sh
# Top replicates together. Can do very stringent filtering for calling peaks 
# 2019-01-15

dmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_K562_round2/all_bam_links"

[[ ! -d $dmain ]] && echo "$dmain not found, exiting" && exit 1

# only two marks are done apparently
marks="H3K27me3 H3K4me1"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged_K562_round2"

[[ ! -d $outdir ]] && mkdir $outdir

for mark in $marks; do
    echo $mark
    inbams=$(echo $dmain/PZ-K562-G1-M-$mark-*.bam)
    outbam=$outdir/K562_"$mark"_merged.bam
    [[ -e $outbam ]] && echo "$outbam exists, skipping" && continue
    samtools merge $outbam $inbams
    samtools index $outbam
done




