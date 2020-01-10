#!/bin/sh
# Jake Yeung
# 1-make_links.sh
# 2019-10-15

inmain="/hpc/hub_oudenaarden/seqdata/OUD3910/AVO516"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3910_HD_0"
[[ ! -d $outmain ]] && mkdir $outmain

for f in `ls -d $inmain/OUD*/PZ*fastq.gz`; do
    fbase=$(basename $f)
    fout=$outmain/$fbase
    ln -s $f $fout
done
