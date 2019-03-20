#!/bin/sh
# Jake Yeung
# 3d-merge_split_files.sh
# After calculating distances, remerge the split files 
# 2019-03-07

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/H3K4me1/motevo_outputs/bed/merged_bed_closestbed_long"
indir=$inmain/split
outf=$inmain/motevo_merged.closest.dist.long.bed

cat $indir/*.bed > $outf
