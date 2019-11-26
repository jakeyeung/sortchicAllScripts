#!/bin/sh
# Jake Yeung
# demultiplex_fastqs.sh
#  
# 2019-09-28

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawDataB6StemCells"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
# grepstr="*.fastq.gz"

cd $inmain
. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3

demux.py *.fastq.gz -hd 0 -merge _ | grep submission | sh

