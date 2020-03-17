#!/bin/sh
# Jake Yeung
# 2019-12-05.copy_B6_other_winsizes.sh
#  
# 2019-12-04

indir="/home/jyeung/hpc/scChiC/raw_data/ZellerRawDataB6_mergedAll.retag/countTables_otherWinSize"
outdir="/home/jyeung/data/from_cluster/scchic/ZellerRawDataB6_mergedAll.retag/countTables_otherWinSize"

scp -r $indir $outdir
