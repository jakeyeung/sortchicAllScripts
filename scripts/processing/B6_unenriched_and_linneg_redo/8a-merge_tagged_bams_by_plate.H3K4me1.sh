#!/bin/sh
# Jake Yeung
# 8a-merge_tagged_bams_by_plate.sh
# Create one tagged bam for each mark 
# 2019-11-26

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_mergedAll"

# Merge H3K4me1: Unenriched and StemCells (lineage neg not yet done)
indir=/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data
infsBM=$indir/ZellerRawDataB6/raw_demultiplexed/PZ-ChIC-Bl6-BM-H3K4me1-Index*-12-09-19/tagged/*.bam
infsSC=$indir/ZellerRawDataB6StemCells/raw_demultiplexed/PZ-ChIC-Bl6-BM-stem-cells-H3K4me1-Index*-12-09-19/tagged/*.bam
outbam=$outdir/H3K4me1-BM_SC-merged.tagged.bam

# [[ ! -e $outbam ]] && echo "$outbam not found, merging" && samtools merge $outbam $infsBM $infsSC --output-fmt BAM
# echo $infsBM $infsSC
samtools merge -f $outbam $infsBM $infsSC --output-fmt BAM
samtools index $outbam


