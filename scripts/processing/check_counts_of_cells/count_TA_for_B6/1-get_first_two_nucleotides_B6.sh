#!/bin/sh
# Jake Yeung
# get_first_two_nucleotides.sh
#  
# 2019-04-16
# must be in t2 

runs="VAN3419 VAN3418 VAN3416"
inmain="/hpc/archive/hub_oudenaarden/seqdata"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_summaries_B6"
[[ ! -d $outdir ]] && mkdir $outdir

n=0
maxjobs=8

# echo "$indir/VAN2979*/PZ-BM-m*H3*R1*.fastq.gz"
for run in $runs; do
    indir=$inmain/$run
    [[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
    for inf in $(ls -d $indir/1904*/Data/Intensities/BaseCalls/AVO*/*/B6-13W1-BM-H3K*.fastq.gz); do
        echo $inf
        # get 4th to 13th basepair
        bname=$(basename $inf)
        bname=${bname%%.*}
        outf=$outdir/$bname.summary.out
        zcat $inf | awk '(NR-2) % 4 == 0' | cut -b 4-13 > $outf&
        if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
            # define maxjobs and n using maxjobsn skeleton
            wait # wait until all have finished (not optimal, but most times good enough)
            echo $n wait
        fi
    done
done

