#!/bin/sh
# Jake Yeung
# 1-merge_bams_from_list.sh
# merge bams from list
# 2019-02-04

jmark="H3K4me1"
bamlistdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bamlists_for_merging/cell_clusters_bin_${jmark}"

wd="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_split_by_bc/count_thres-0"

cd $wd

n=0
maxjobs=4

for blist in `ls -d $bamlistdir/bamlist.*.txt`; do
    # echo $blist
    clstr=`echo $blist | cut -d"." -f2`
    # echo $clstr
    bout="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_cluster_merged/${jmark}_cluster_${clstr}.bam"
    echo "samtools merge -b $blist $bout"
    samtools merge -b $blist $bout&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait
