#!/bin/sh
# Jake Yeung
# 3-download_raw_Bartosovic.sh
#  
# 2021-06-07

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data"
inf="${indir}/metadata_Wu_noheader.txt"
outdir=${indir}/prefetch_outputs

[[ ! -d $outdir ]] && mkdir $outdir

srx="SRR12638101"
prefetch -v $srx -O $outdir
