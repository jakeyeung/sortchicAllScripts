#!/bin/sh
# Jake Yeung
# 2-merge_beds_together_CUTnTAG.sh
#  
# 2021-06-05
# Merge beds together to do peak calling 

jmark="H3K4me2"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_Kaya-Okur/${jmark}"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_Kaya-Okur/${jmark}/merged_beds"
[[ ! -d $outdir ]] && mkdir $outdir

outbed="${outdir}/K562_${jmark}_20181120_allmerged.bed"
zcat ${indir}/*.bed.gz > $outbed
