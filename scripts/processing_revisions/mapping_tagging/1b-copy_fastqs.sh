#!/bin/sh
# Jake Yeung
# 1b-copy_fastqs.sh
#  
# 2022-01-21

indir="/hpc/archive/hub_oudenaarden/seqdata/OUD6998"
indir2="/hpc/archive/hub_oudenaarden/seqdata/OUD6999"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/revisions_data/new_experiments/raw_data"

scp -r $indir $indir2 $outdir
