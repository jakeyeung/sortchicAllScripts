#!/bin/sh
# Jake Yeung
# 3c-filter_bed_by_distance.sh
# Filter by distance 
# 2019-03-07

# WRAP UP
while [[ `qstat | grep "Annot_" | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

threses="0 1000 10000"
marks="H3K4me1"

n=0
maxjobs=3
for thres in $threses; do
    for mark in $marks; do
        inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_cluster_BM-AllMerged_Peaks_1000/${jmark}/motevo_outputs/bed/merged_bed_closestbed_long"
        bname="motevo_merged.closest.long"
        inf=${inmain}/${bname}.reannotated.bed
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
