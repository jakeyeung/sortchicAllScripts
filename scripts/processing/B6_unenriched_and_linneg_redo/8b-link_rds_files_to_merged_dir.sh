#!/bin/sh
# Jake Yeung
# 8b-link_rds_files_to_merged_dir.sh
#  
# 2019-11-23

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6_all/countTables"
[[ ! -d $outdir ]] && mkdir $outdir

for f in `ls -d $inmain/ZellerRawDataB6*/countTables/*demuxbugfixed.rds`; do
    # ffull=$(readlink -e $f)
    # echo $ffull
    ln -s $f $outdir/.
done
