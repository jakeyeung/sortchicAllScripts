#!/bin/sh
# Jake Yeung
# 3b-reannotate_distance_on_bed.sh
# Reannotate distance on bed 
# 2019-03-07

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/annotate_distance_on_bed.R"
wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
cd $wd

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/$jmark/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.long.bed"
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/H3K4me1/motevo_outputs/bed/merged_bed_closestbed_long/test.bed"
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/H3K4me1/motevo_outputs/bed/merged_bed_closestbed_long/test.dist.bed"
# outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/$jmark/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.dist.long.bed"

[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
[[ -e $outf ]] && echo "$outf found, exiting" && exit 1

Rscript $rs $inf $outf

# for jmark in $jmarks; do
# 
#     inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/$jmark/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.long.bed"
#     outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/$jmark/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.dist.long.bed"
# 
#     [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
#     [[ -e $outf ]] && echo "$outf found, exiting" && exit 1
# 
#     Rscript $rs $inf $outf&
# done
# wait
