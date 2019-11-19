#!/bin/sh
# Jake Yeung
# 1-demultiplex.sh
#  
# 2019-09-04

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/fastq"
grepstr="PZ-ChIC*.fastq.gz"

cd $inmain
. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
echo "demux.py "${grepstr}" -merge _ | grep submission | sh"
