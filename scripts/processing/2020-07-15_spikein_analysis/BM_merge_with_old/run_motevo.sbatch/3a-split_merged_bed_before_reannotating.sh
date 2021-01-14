#!/bin/sh
# Jake Yeung
# 3b-split_merged_bed_before_reannotating.sh
# Annotate distances again but merged file too large. Split to smaller files 
# 2019-03-07

# mark="H3K9me3"
# mark="H3K4me1"
mark="H3K4me3"
# mark="H3K27me3"

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_ZF-AllMerged2_Peaks_1000"
inf="${inmain}/${mark}/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.long.bed"
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

outmain=$(dirname $inf)
outdir=$outmain/split
[[ ! -d $outdir ]] && mkdir $outdir
split -l 10000000 $inf $outdir/motevo_merged_split.

