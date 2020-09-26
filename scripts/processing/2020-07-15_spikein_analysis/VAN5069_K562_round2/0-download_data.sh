#!/bin/sh
# Jake Yeung
# 0-download_data.sh
#  
# 2020-09-07

link="https://ncie01.op.umcutrecht.nl/s/XswAaie8erZRsXS/download"
# link="https://ncie01.op.umcutrecht.nl/s/XswAaie8erZRsXS/download"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2"

cd $indir

curl -O $link

# rename to VAN5039.tar.gpg manually 
