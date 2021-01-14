#!/bin/sh
# Jake Yeung
# 2b-unzip_split.sh
#  
# 2020-02-16

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM-AllMerged_Peaks_1000/H3K4me1/motevo_outputs"
cd $indir
# strip 8 because we want to start at split 
# hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM-AllMerged_Peaks_1000/H3K4me1/motevo_outputs/split/ar/Ebf3.pwm.prior.gz

tar -zxvf split.tar.gz --strip=8
