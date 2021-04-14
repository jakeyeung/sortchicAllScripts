#!/bin/sh
# Jake Yeung
# 0-copy_fastq_to_dir.sh
#  
# 2020-10-03

inmain="/hpc/archive/hub_oudenaarden/seqdata/VAN5109"
cd $inmain
find . -name "*PZ*.fastq.gz" -exec cp {} /hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5109_BM_twoplates/. \;

