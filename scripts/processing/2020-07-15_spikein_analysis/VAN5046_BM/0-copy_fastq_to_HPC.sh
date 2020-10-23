#!/bin/sh
# Jake Yeung
# 0-copy_fastq_to_HPC.sh
#  
# 2020-09-11

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5046_BM"
link="https://ncie01.op.umcutrecht.nl/s/NqLscSBitw5yD6n/download"

cd $indir 

curl -O $link
