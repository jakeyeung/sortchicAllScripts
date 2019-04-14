#!/bin/sh
# Jake Yeung
# 4-bam_to_bigwig.sh
#  
# 2019-02-04

jdate="2019-03-28"
suffix="build95"

bs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/lib/bam_to_bigwig_mm10.sh"
[[ ! -e $bs ]] && echo "$bs not found, exiting" && exit 1
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/sorted_bams_${suffix}_${jdate}"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bigwigs_${suffix}_${jdate}"

[[ ! -d $outdir ]] && mkdir $outdir

bsize=1000

n=0
maxjobs=16
for b in `ls -d $indir/*.bam`; do
        bbase=$(basename $b)
        bbase=${bbase%%.*}
        # echo $bbase
        bout=$outdir/$bbase.bw
        # echo "bash $bs $b $bout&"
        bash $bs $b $bout $bsize&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait
