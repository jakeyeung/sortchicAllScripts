#!/bin/sh
# Jake Yeung
# 2019-03-21_send_public_bigwigs.sh
#  
# 2019-03-21

indir1="/Users/yeung/data/scchic/public_data/Bcells"
indir2="/Users/yeung/data/scchic/public_data/Neutrophils"
outdir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data"

scp -r $indir1 $indir2 $outdir
