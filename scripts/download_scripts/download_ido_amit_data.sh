#!/bin/sh
# Jake Yeung
# download_ido_amit_data.sh
# Download ido amit data iCHIP 
# 2019-03-25

inf="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60103/suppl//GSE60103_RAW.tar"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Lara-Astiaso_2014_Science"

cd $indir

curl -O $inf
