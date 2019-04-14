#!/bin/sh
# Jake Yeung
# run.variance_over_pseudotime_plots_parallel_for_server.sh
#  
# 2019-04-04

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/scripts_analysis/build95_primetime_from_server/variance_over_pseudotime_plots_parallel_for_server.R"
wd="/home/hub_oudenaarden/jyeung/projects/scChiC"

cd $wd; Rscript $rs 
