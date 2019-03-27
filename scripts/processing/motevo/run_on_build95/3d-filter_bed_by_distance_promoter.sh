#!/bin/sh
# Jake Yeung
# 3c-filter_bed_by_distance.sh
# Filter by distance 
# 2019-03-07

thres=1000

jmark="H3K4me1"
# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/H3K4me3/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.dist.long.bed"
# outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/H3K4me3/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.dist.filt.$thres.promoter.long.bed"
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/$jmark/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.dist.long.bed"
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/$jmark/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.dist.filt.$thres.promoter.long.bed"


awk '$7 < 1000 && $7 > -1000' $inf > $outf
