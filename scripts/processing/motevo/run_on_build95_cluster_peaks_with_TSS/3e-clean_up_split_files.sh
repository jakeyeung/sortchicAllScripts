#!/bin/sh
# Jake Yeung
# 3e-clean_up_split_files.sh
# Clean up splits after reannotating promtoer  
# 2019-03-26

# mark="H3K27me3"
mark="H3K4me3"
suffix=""
indir1="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene_cluster_build95/$mark/motevo_outputs/bed/merged_bed_closestbed_long/split"
# indir2="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene_build95/H3K4me3/motevo_outputs/bed/merged_bed_closestbed_long/split"

[[ ! -d $indir1 ]] && echo "$indir1 not found, exiting" && exit 1

echo "Run these manually for safety"
echo "rm $indir1"
