#!/bin/sh
# Jake Yeung
# 2019-12-05.copy_giladi_markergenes.sh
#  
# 2019-12-06

inf="$HOME/hpc/scChiC/public_data/Giladi_et_al_2018/41556_2018_121_MOESM4_ESM.markergenes.xlsx"
outf="$HOME/data/from_cluster/public_data/Giladi_et_al_2018/."

cp $inf $outf
