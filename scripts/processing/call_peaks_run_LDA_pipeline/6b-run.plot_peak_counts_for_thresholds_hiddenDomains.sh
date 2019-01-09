#!/bin/sh
# Jake Yeung
# 6b-run.plot_peak_counts_for_thresholds.sh
# Run plot peak counts 
# 2018-12-22

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/plot_peak_counts_for_thresholds.R"

# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats.merge_25000/PZ-BM-H3K4me1.merged.NoCountThres.merge_25000.Robj"  # output of make count mat


marks="H3K27me3 H3K4me1 H3K4me3 H3K9me3"
for mark in $marks; do
    echo $mark
    inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromHiddenDomains.1000/PZ-BM-${mark}.merged.NoCountThres.hiddenDomains.Robj"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats.fromHiddenDomains.1000/threshold_plots.${mark}"
    [[ ! -d $outdir ]] && mkdir $outdir
    cd $wd; Rscript $rs $inf $outdir
done


