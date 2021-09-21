#!/bin/sh
# Jake Yeung
# 1-decypt_gpg.sh
#  
# 2020-09-07

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5232_VAN5233_BM"

cd $indir
# gpg VAN5039.tar.gpg  /  vanoudenaarden
gpg VAN5233.tar.gpg
