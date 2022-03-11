#!/bin/sh
# Jake Yeung
# 2-merge_beds_together_CUTnTAG.sh
#  
# 2021-06-05
# Merge beds together to do peak calling 


indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_Kaya-Okur/H3K27me3"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_Kaya-Okur/H3K27me3/merged_beds"

outbed="${outdir}/K562_H3K27me3_20181120_allmerged.bed"
zcat ${indir}/*.bed.gz > $outbed
