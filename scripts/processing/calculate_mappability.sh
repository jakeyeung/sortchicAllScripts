#!/bin/sh
# Jake Yeung
# calculate_mappability.sh
#  
# 2019-06-10

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
reps="1 2 3 4"

dirmain="/hpc/hub_oudenaarden/avo/scChiC/raw_demultiplexed"

for mark in $marks; do
    for rep in $reps; do
        indir=$dirmain/"B6-13W1-BM-$mark-$rep-merged"
        [[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
        bamin="$indir/bwaMapped.bam"
        [[ ! -e $bamin ]] && echo "$bamin not found, exiting" && exit 1
        echo "samtools flagstat $bamin"
        samtools flagstat $bamin
    done
done
