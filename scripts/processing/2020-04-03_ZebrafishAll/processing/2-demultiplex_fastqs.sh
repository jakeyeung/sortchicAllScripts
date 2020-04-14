#!/bin/sh
# Jake Yeung
# demultiplex_fastqs.sh
#  
# 2019-09-28

inmain="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/raw_data"

cd $inmain
. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3

demux.py *.fastq.gz -hd 0 | grep submission | sh
