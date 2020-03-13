#!/bin/sh
# Jake Yeung
# copy_counts_and_RZ_BM_all_demuxbugfix_retagged.sh
# Copy things to Rstudio 
# 2019-12-03

rzdir="scChiC/raw_data/ZellerRawDataB6_mergedAll.retag/RZcounts"
countdir="scChiC/raw_data/ZellerRawDataB6_mergedAll.retag/countTables"
outdir="/home/jyeung/data/from_cluster/scchic/ZellerRawDataB6_mergedAll.retag"

mkdir $outdir

cp -r ~/hpc/$rzdir ~/hpc/$countdir $outdir
