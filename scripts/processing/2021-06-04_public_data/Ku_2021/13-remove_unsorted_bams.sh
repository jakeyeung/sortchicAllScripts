#!/bin/sh
# Jake Yeung
# 13-remove_unsorted_bams.sh
#  
# 2021-06-15

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Ku_et_al_2021/SRA_data/prefetch_outputs/bams_demux_bugfix"
cd $indir
cmd="rm *_1.bam"
echo $cmd
