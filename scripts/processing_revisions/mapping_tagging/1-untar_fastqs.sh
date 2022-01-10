#!/bin/sh
# Jake Yeung
# 1-untar_fastqs.sh
#  
# 2022-01-06

indir="/hpc/hub_oudenaarden/Peter/seqdata"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/fastqs"

cd $indir
for f in PZ*.tar; do
    echo $f
    tar -xvf $f --directory $outdir 
done
