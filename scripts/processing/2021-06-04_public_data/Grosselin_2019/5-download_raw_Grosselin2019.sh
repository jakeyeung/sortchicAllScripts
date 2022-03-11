#!/bin/sh
# Jake Yeung
# 3-download_raw_Bartosovic.sh
#  
# 2021-06-07

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/SRA_data"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Grosselin_et_al_2019/SRA_data"
inf="${indir}/metadata_Grosselin2019_noheader.txt"
outdir=${indir}/prefetch_outputs

[[ ! -d $outdir ]] && mkdir $outdir

cd $outdir

while read p; do
  srx=`echo "$p" | cut -f19`
  echo $srx
  prefetch -v $srx -O $outdir
  # exit 0
done < $inf
