#!/bin/sh
# Jake Yeung
# 5-run.select_best_gene.sh
#  
# 2019-04-19

n=0
maxjobs=4

suffix="0_build95_B6.withchr"
tssdist=50000
marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
prefix="B6-BM"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/call_peaks_on_bam_clusters/select_best_gene.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

for mark in $marks; do

    inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromGeneTSS.${suffix}_${tssdist}.cells_from_bin_analysis/${prefix}-${mark}.merged.NoCountThres.GeneTSS.Robj"
    outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromGeneTSS.${suffix}_${tssdist}.cells_from_bin_analysis/${prefix}-${mark}.merged.NoCountThres.GeneTSS.Dedup.Robj"

    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    [[ -e $outf ]] && echo "$outf found, exiting" && exit 1
    Rscript $rs $inf $outf&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
    	# define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait

