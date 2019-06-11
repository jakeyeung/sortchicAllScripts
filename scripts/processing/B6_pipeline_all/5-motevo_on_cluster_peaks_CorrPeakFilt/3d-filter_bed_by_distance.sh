#!/bin/sh
# Jake Yeung
# 3c-filter_bed_by_distance.sh
# Filter by distance 
# 2019-03-07

prefix="cluster"
suffix="build95_B6_CorrPeakFilt"
threses="0 1000 10000"
# marks="H3K4me1 H3K4me3"
# marks="H3K27me3 H3K9me3"
# marks="H3K4me1"
# marks="H3K4me3"
marks="H3K27me3"
# marks="H3K27me3"
# marks="H3K9me3"

n=0
maxjobs=3
for thres in $threses; do
    for mark in $marks; do
        inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene_${prefix}_${suffix}/${mark}/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.dist.long.bed"
        outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene_${prefix}_${suffix}/${mark}/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.dist.filt.${thres}.long.bed"
        [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
        awk -v jthres=$thres '$7 > jthres || $7 < -jthres' $inf > $outf&
        if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
            # define maxjobs and n using maxjobsn skeleton
            wait # wait until all have finished (not optimal, but most times good enough)
            echo $n wait
        fi
    done

done
wait
