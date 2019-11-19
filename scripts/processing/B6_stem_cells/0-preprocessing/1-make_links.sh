#!/bin/sh
# Jake Yeung
# 1-make_links.sh
#  
# 2019-09-28

inmain="/hpc/hub_oudenaarden/jyeung/raw_data_from_sequencer/AVO508"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3700"
[[ ! -d $outmain ]] && mkdir $outmain

for f in `ls -d $inmain/OUD*/*fastq.gz`; do
    fbase=$(basename $f)
    fout=$outmain/$fbase
    ln -s $f $fout
done
