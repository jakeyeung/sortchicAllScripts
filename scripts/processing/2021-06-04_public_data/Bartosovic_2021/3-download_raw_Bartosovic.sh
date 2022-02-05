#!/bin/sh
# Jake Yeung
# 3-download_raw_Bartosovic.sh
#  
# 2021-06-07

# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bartosovic_et_al_2021/SRA_data/sra_result_Bartosovic_noheader.csv"
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bartosovic_et_al_2021/SRA_data/metadata_Bartosovic_noheader.txt"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bartosovic_et_al_2021/SRA_data/prefetch_outputs"
[[ ! -d $outdir ]] && mkdir $outdir

while read p; do
  srx=`echo "$p" | cut -f19`
  echo $srx
  prefetch -v $srx -O $outdir
  # exit 0
done < $inf
