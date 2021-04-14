#!/bin/sh
# Jake Yeung
# 0-download_data.sh
#  
# 2020-10-05

# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5232_VAN5233_BM"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5234_VAN5235_VAN5236_BM"
cd $outdir

outname="VAN5234.tar.gpg"

# link="https://ncie01.op.umcutrecht.nl/s/7bqWtCSm2GWxn8E/download"
link="https://ncie01.op.umcutrecht.nl/s/WmXnJF56mzaTmQa/download"

wget --output-document=${outname} $link
