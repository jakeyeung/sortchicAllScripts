#!/bin/sh
# Jake Yeung
# 3b-split_merged_bed_before_reannotating.sh
# Annotate distances again but merged file too large. Split to smaller files 
# 2019-03-07


# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/H3K4me1/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.long.bed"
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/H3K27me3/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.long.bed"

outmain=$(dirname $inf)
outdir=$outmain/split
[[ ! -d $outdir ]] && mkdir $outdir
split -l 10000000 $inf $outdir/motevo_merged_split.

