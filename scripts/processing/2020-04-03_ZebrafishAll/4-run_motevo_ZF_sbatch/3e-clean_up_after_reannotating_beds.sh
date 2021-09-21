#!/bin/sh
# Jake Yeung
# 3e-clean_up_after_reannotating_beds.sh
#  
# 2020-08-25

jmark="H3K4me1"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_ZF-AllMerged2_Peaks_1000/${jmark}/motevo_outputs/bed/merged_bed_closestbed_long/split"

cmd="Remove this dir: $indir"
echo $cmd
