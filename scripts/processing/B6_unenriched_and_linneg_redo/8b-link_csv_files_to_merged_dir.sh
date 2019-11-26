#!/bin/sh
# Jake Yeung
# 8b-link_csv_files_to_merged_dir.sh
#  
# 2019-11-22

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_all/LHcounts"
[[ ! -d $outdir ]] && mkdir $outdir

for f in `ls -d $inmain/ZellerRawDataB6*/RZcounts/*LH_counts.demuxbugfixed.csv`; do
    # ffull=$(readlink -e $f)
    # echo $ffull
    ln -s $f $outdir/.
done
