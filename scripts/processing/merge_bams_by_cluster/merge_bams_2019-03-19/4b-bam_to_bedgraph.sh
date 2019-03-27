#!/bin/sh
# Jake Yeung
# 4-bam_to_bigwig.sh
#  
# 2019-02-04

bs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/lib/bam_to_bedgraph_mm10.sh"

[[ ! -e $bs ]] && echo "$bs not found, exiting" && exit 1
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/sorted_bams_2019-03-20"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/bedgraphs_2019-03-21"
[[ ! -d $outdir ]] && mkdir $outdir

n=0
maxjobs=16
for b in `ls -d $indir/*.bam`; do
        bbase=$(basename $b)
        bbase=${bbase%%.*}
        # echo $bbase
        bout=$outdir/$bbase.bedgraph
        bash $bs $b $bout&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait
