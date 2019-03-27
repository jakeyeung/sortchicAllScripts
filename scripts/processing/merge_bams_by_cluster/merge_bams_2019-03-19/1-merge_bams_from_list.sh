#!/bin/sh
# Jake Yeung
# 1-merge_bams_from_list.sh
# merge bams from list
# 2019-03-19

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/sorted_bams_2019-03-20"
[[ ! -d $outdir ]] && mkdir $outdir
wd="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_split_by_bc/count_thres-0"
[[ ! -d $wd ]] && echo "$wd not found, exiting" && exit 1
cd $wd

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

bamlistdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/bamlist_for_merging"

n=0
maxjobs=8

    for blist in `ls -d $bamlistdir/bamlist-*.txt`; do
        bname=$(basename $blist)
        bname=${blist%%.*}
        clstr=`echo $bname | cut -d"-" -f2`
        jmark=`echo $bname | cut -d"-" -f3`
        bout="$outdir/${jmark}_cluster_${clstr}.bam"
        # echo "samtools merge -b $blist $bout"
        samtools merge -b $blist $bout&
        if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
            # define maxjobs and n using maxjobsn skeleton
            wait # wait until all have finished (not optimal, but most times good enough)
            echo $n wait
        fi
    done
wait
