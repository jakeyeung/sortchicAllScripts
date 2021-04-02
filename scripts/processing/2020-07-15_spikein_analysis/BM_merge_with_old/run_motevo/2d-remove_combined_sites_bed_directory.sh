#!/bin/sh
# Jake Yeung
# 2d-remove_combined_sites_bed_directory.sh
#  
# 2020-11-10

jmarks="H3K4me1 H3K4me3 H3K27me3"

for jmark in $jmarks; do
    indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM-AllMerged3_Peaks/${jmark}/motevo_outputs/bed"
    cd $indir
    dname="combined_sites_bed"
    cmd="rm -r ${indir}/${dname}"
    echo $cmd
done
