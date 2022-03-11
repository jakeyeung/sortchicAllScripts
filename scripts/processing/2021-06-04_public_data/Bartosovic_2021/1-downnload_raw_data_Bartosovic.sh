#!/bin/sh
# Jake Yeung
# 1-downnload_raw_data_Bartosovic.sh
#  
# 2021-06-07

# download everhything

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bartosovic_et_al_2021"
cd $indir
link="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163532/suppl/GSE163532_RAW.tar"

curl -O $link


