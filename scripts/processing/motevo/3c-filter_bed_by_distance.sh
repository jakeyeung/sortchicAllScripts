#!/bin/sh
# Jake Yeung
# 3c-filter_bed_by_distance.sh
# Filter by distance 
# 2019-03-07

thres=10000

mark="H3K4me1"
# mark="H3K4me3"

inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/$mark/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.dist.long.bed"
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/$mark/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.dist.filt.$thres.long.bed"
# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/H3K4me1/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.dist.long.bed"
# outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/H3K4me1/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.dist.filt.$thres.long.bed"

[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

awk -v jthres=$thres '$7 > jthres || $7 < -jthres' $inf > $outf
