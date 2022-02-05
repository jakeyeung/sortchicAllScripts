#!/bin/sh
# Jake Yeung
# download_data_from_GEO.sh
#  
# 2021-06-07

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021"
link="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE139nnn/GSE139857/suppl"

wget -nd -r -P $outdir $link
