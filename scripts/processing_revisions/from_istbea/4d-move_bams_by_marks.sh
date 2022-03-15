#!/bin/sh
# Jake Yeung
# 4d-move_bams_by_marks.sh
#  
# 2022-01-25

jmarks="k4me1 k4me3 k27me3 k9me3"
inbase="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/tagged_bams"
cd $inbase

for jmark in $jmarks; do
  outdir="${inbase}/BM_${jmark}"
  [[ ! -d $outdir ]] && mkdir $outdir
  # echo "mv *${jmark}* ${outdir}"
  # ls *${jmark}*
  mv PZ*${jmark}* ${outdir}
done
