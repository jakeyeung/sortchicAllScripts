#!/bin/sh
# Jake Yeung
# 2019-01-22_copy_binned_matrix.sh
#  
# 2019-01-22

inftxt="/hpc/hub_oudenaarden/avo/scChiC/metacell/BM-*-100kb.txt"
outdir="/tmp/binned_mat/."

scp t2:$inftxt $outdir
