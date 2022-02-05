#!/bin/sh
# Jake Yeung
# 5-prefetch_from_metadata.sh
#  
# 2021-06-09

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data"
inf="${indir}/metadata_Ku2021.noheader.txt"
outdir=${indir}/prefetch_outputs

[[ ! -d $outdir ]] && mkdir $outdir
cd $outdir

while read p; do
  srx=`echo "$p" | cut -f19`
  echo $srx
  prefetch -v $srx -O $outdir
  # exit 0
done < $inf
