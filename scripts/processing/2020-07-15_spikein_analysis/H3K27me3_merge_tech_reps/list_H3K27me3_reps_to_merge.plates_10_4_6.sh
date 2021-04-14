#!/bin/sh
# Jake Yeung
# list_H3K27me3_reps_to_merge.sh plates 10, 4, 6 are in VAN5046 without rep3 name
#  
# 2020-11-21

indir="/hpc/archive/hub_oudenaarden/seqdata/VAN5046/200910_NS500813_0646_AHFLJGBGXG/Data/Intensities/BaseCalls/AVOEI856-27"
cd $indir

find . -name "PZ-BM*rep3-H3K27me3-*fastq.gz" 2>/dev/null


