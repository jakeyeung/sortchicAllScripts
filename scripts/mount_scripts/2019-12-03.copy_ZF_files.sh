#!/bin/sh
# Jake Yeung
# 2019-12-03.copy_ZF_files.sh
#  
# 2019-12-03

indir1="/home/jyeung/hpc/scChiC/raw_data/ZellerRawDataZF_all.retag/countTables"
indir2="/home/jyeung/hpc/scChiC/raw_data/ZellerRawDataZF_all.retag/RZcounts"
outdir="/home/jyeung/data/from_cluster/scchic/ZellerRawDataZF_all.retag"

cp -r $indir1 $indir2 $outdir


