#!/bin/sh
# Jake Yeung
# demultiplex_fastqs.sh
#  
# 2019-09-28

# inmain="/hpc/hub_oudenaarden/jyeung/raw_data_from_sequencer/AVO508/links_merged"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData"
# grepstr="PZ-ChIC*.fastq.gz"
grepstr="*fastq.gz"

cd $inmain
. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
demux.py "${grepstr}" -merge _ | grep submission | sh

