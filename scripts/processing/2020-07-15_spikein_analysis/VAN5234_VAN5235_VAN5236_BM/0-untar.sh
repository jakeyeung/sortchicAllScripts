#!/bin/sh
# Jake Yeung
# 0-untar.sh
#  
# 2020-10-05

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5232_VAN5233_BM"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5234_VAN5235_VAN5236_BM"
cd $indir

tar -xvf "VAN5234.tar"&
tar -xvf "VAN5235.tar"&
tar -xvf "VAN5236.tar"&
wait
