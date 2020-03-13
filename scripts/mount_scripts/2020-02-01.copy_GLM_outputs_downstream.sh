#!/bin/sh
# Jake Yeung
# 2020-02-01.copy_GLM_outputs.sh
#  
# 2020-02-02

indir="/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.downstream/"
outdir="/home/jyeung/hpc/scChiC/from_rstudioserver/GLMPCA_outputs.downstream"

rsync -avrL --copy-links $indir $outdir
# scp $indir/*  $outdir

