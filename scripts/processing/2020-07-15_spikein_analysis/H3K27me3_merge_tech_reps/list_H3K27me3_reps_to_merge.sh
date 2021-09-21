#!/bin/sh
# Jake Yeung
# list_H3K27me3_reps_to_merge.sh
#  
# 2020-11-21

indir="/hpc/archive/hub_oudenaarden/seqdata"
cd $indir

find . -name "PZ-BM*rep3-H3K27me3-*fastq.gz" 2>/dev/null


