#!/bin/sh
# Jake Yeung
# 1-copy_OUD3913_H3K27me3_to_OUD3985.sh
#  
# 2019-11-12

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/fastq_tmp"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3985"

for f in `ls -d $inmain/PZ*.fastq.gz`; do
    ln -s $f $outdir
done
