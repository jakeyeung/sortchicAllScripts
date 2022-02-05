#!/bin/sh
# Jake Yeung
# download_data_from_GEO.sh
#  
# 2021-06-07

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021"
link="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157910/suppl/GSE157910_RAW.tar"

cd $indir
curl -O $link
