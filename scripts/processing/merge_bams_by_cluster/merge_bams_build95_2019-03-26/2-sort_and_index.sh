#!/bin/sh
# Jake Yeung
# 2-sort_and_index.sh
# After merging, sort and index 
# 2019-02-04

suffix="build95"
jdate="2019-03-26"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/sorted_bams_${suffix}_${jdate}"

n=0
maxjobs=4

for bsorted in `ls -d $indir/*.bam`; do
    samtools index $bsorted&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait
