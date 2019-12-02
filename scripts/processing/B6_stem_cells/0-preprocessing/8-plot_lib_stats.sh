#!/bin/sh
# Jake Yeung
# 7-plot_lib_stats.sh
#  
# 2019-09-04

inmain="/hpc/hub_oudenaarden/jyeung/raw_data_from_sequencer/AVO508/links_merged/raw_demultiplexed"
cd $inmain
libs=$(for d in `ls -d $inmain/PZ-ChIC*`; do echo $(basename $d); done)

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3; libraryStatistics.py $libs  -tagged_bam /tagged/bwaMapped.sorted.bam

