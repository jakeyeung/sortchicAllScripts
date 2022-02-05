#!/bin/sh
# Jake Yeung
# download_data_from_GEO.sh
#  
# 2021-06-07

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Grosselin_et_al_2019"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2019"
[[ ! -d $outdir ]] && mkdir $outdir
# link="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE105nnn/GSE105012/suppl"
# GSE117309 
link="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117309/suppl/GSE117309_RAW.tar"
cd $outdir
curl -O $link
# wget -nd -r -P $outdir $link
