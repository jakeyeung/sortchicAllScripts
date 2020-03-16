#!/bin/sh
# Jake Yeung
# 2019-12-05.copy_giladi.sh
#  
# 2019-12-05

indir="$HOME/hpc/scChiC/public_data/Giladi_et_al_2018"
outdir="/home/jyeung/data/from_cluster/public_data"
[[ ! -d $outdir ]] && mkdir $outdir

cp -r $indir $outdir
