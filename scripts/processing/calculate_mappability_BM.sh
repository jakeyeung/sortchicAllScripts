#!/bin/sh
# Jake Yeung
# calculate_mappability.sh
#  
# 2019-06-10

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
reps="m1 m2"
techreps="1 2"
# reps="rep1_H2VTWBGX9_S rep2_H2VTWBGX9_S rep3_AH2GNLBGX9_S"
# reps=

dirmain="/hpc/hub_oudenaarden/avo/scChiC/raw_demultiplexed"

for mark in $marks; do
    for rep in $reps; do
        for techrep in $techreps; do
        indirprefix=$dirmain/"PZ-BM-$rep-$mark-$techrep"
            for indir in $indirprefix*; do
                [[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
                bamin="$indir/bwaMapped.bam"
                [[ ! -e $bamin ]] && echo "$bamin not found, exiting" && exit 1
                echo "samtools flagstat $bamin"
                samtools flagstat $bamin
            done
        done
    done
done
