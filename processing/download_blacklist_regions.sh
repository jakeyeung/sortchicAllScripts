#!/bin/sh
# Jake Yeung
# download_blacklist_regions.sh
# Download blacklist regions 
# 2018-12-18

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/blacklist"

cd $outdir

curl -O http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz 
