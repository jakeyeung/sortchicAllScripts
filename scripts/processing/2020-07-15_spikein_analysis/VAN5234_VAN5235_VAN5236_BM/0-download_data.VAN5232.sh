#!/bin/sh
# Jake Yeung
# 0-download_data.sh
#  
# 2020-10-05

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5232_VAN5233_BM"
cd $outdir

outname="VAN5232.tar.gpg"

[[ -e $outdir/$outname ]] && echo "$outdir/$outname found, exiting" && exit 1

# 
# link="https://ncie01.op.umcutrecht.nl/s/7bqWtCSm2GWxn8E/download"
link="https://ncie01.op.umcutrecht.nl/s/L4LFS8NYJ26WWTi/download"

wget --output-document=${outname} $link
