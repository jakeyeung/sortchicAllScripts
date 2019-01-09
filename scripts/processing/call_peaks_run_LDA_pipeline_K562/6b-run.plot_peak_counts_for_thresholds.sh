#!/bin/sh
# Jake Yeung
# 6b-run.plot_peak_counts_for_thresholds.sh
# Run plot peak counts 
# 2018-12-22

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/plot_peak_counts_for_thresholds.R"

inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats.merge_25000/PZ-BM-H3K4me1.merged.NoCountThres.merge_25000.Robj"  # output of make count mat
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats.merge_25000/plots"

[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
[[ ! -d $outdir ]] && mkdir $outdir

cd $wd; Rscript $rs $inf $outdir
