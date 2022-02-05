#!/bin/sh
# Jake Yeung
# download_brains.sh
#  
# 2021-06-07

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bartosovic_et_al_2021"
link="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE157nnn/GSE157637/suppl"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Bartosovic_et_al_2021/GSE157637"

cd $indir 

# curl -O $link
wget -nd -r -P $outdir $link
