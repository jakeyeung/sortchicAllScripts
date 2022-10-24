#!/bin/sh
# Jake Yeung
# 2-index_bams.sh
#  
# 2022-07-18

indir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/bams"

n=0
maxjobs=6

module load samtools
for f in `ls -d $indir/*.bam`; do
  echo $f
  samtools index $f&
  if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
  	# define maxjobs and n using maxjobsn skeleton
      wait # wait until all have finished (not optimal, but most times good enough)
      echo $n wait
  fi
done
wait
