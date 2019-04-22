#!/bin/sh
# Jake Yeung
# get_first_two_nucleotides.sh
#  
# 2019-04-16

# must be in t2 

indir="/hpc/archive/hub_oudenaarden/seqdata/VAN2979/181107_NS500413_0510_AH3VGVBGX9/Data/Intensities/BaseCalls/AVO443"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_summaries_VAN2979"
[[ ! -d $outdir ]] && mkdir $outdir

n=0
maxjobs=8

# echo "$indir/VAN2979*/PZ-BM-m*H3*R1*.fastq.gz"
for inf in $(ls -d $indir/VAN2979*/PZ-BM-m*H3*R1*.fastq.gz); do
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

