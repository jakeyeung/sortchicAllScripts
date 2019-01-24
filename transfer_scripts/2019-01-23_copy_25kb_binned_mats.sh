#!/bin/sh
# Jake Yeung
# 2019-01-23_copy_25kb_binned_mats.sh
# Copy data for 25kb bins 
# 2019-01-23

indir="/Users/yeung/Dropbox/scCHiC_figs/25kb_bins"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_dropbox/."

scp -r $indir t2:$outdir
