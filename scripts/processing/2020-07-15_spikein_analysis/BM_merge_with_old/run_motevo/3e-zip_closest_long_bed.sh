#!/bin/sh
# Jake Yeung
# 3e-zip_closest_long_bed.sh
#  
# 2020-11-10

jmarks="H3K4me1 H3K4me3 H3K27me3"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM-AllMerged3_Peaks"
suffix="motevo_outputs/bed/merged_bed_closestbed_long"

for jmark in $jmarks; do
    indir2=$indir/$jmark/${suffix}
    cd $indir2
    fname="motevo_merged.closest.long.bed"
    [[ ! -e ${fname} ]] && echo "${fname} not found, exiting" && exit 1
    echo $indir2/${fname}
    gzip $fname
done
