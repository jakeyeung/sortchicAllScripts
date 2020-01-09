#!/bin/sh
# Jake Yeung
# 8a-merge_tagged_bams_by_plate.sh
# Create one tagged bam for each mark 
# 2019-11-26

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_mergedAll"
indir=/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

jmark="H3K9me3"
infsBM=$indir/ZellerRawDataB6/raw_demultiplexed/B6-13W1-BM-${jmark}-*/tagged/*.bam
infsLN=$indir/ZellerRawDataB6LinNeg/raw_demultiplexed/PZ-Bl6-BM-Linneg-${jmark}-*/tagged/*.bam
# infsSC=$indir/ZellerRawDataB6StemCells/raw_demultiplexed/PZ-ChIC-B6BMSC-${jmark}-*/tagged/*.bam  # non yet?? still sequencing??
outbam=$outdir/${jmark}-BM_SC-merged.tagged.bam
[[ ! -e $outbam ]] && echo "$outbam not found, merging" && samtools merge $outbam $infsBM $infsLN --output-fmt BAM


