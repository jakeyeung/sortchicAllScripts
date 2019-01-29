#!/bin/sh
# Jake Yeung
# 2019-01-15_copy_count_mat_K562.sh
#  
# 2019-01-15

# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_mat_K562/count_mats.hiddenDomains.1000"

inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_mat_K562/count_mats.hiddenDomains.1000_round2"

outdir="/tmp/count_mat_K562_round2"
[[ ! -d $outdir ]] && mkdir $outdir

cd $outdir

scp -r t2:$inf $outdir
