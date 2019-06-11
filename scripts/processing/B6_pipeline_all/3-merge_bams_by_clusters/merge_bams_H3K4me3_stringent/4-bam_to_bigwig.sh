#!/bin/sh
# Jake Yeung
# 4-bam_to_bigwig.sh
#  
# 2019-02-04

# jdate="2019-05-09"
suffix="build95_B6_stringent"

bs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/lib/bam_to_bigwig_mm10.sh"
[[ ! -e $bs ]] && echo "$bs not found, exiting" && exit 1
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/sorted_bams_${suffix}"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bigwigs_${suffix}"

[[ ! -d $outdir ]] && mkdir $outdir
bsize=1000
n=0
maxjobs=8
for b in `ls -d $indir/*.bam`; do
        bbase=$(basename $b)
        bbase=${bbase%%.*}
        # echo $bbase
        bout=$outdir/$bbase.bw
        [[ -e $bout ]] && echo "$bout not found, continuing" && continue
        # echo "bash $bs $b $bout&"
        bash $bs $b $bout $bsize&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait
