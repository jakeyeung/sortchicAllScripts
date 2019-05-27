#!/bin/sh
# Jake Yeung
# run.make_fig4_chromosome_over_traj_server_paths.sh
#  
# 2019-05-27

rs="/hpc/hub_oudenaarden/jyeung/code_for_analysis/scchic/scripts/scripts_analysis/B6_analysis/make_figs/make_fig4_chromosome_over_traj_server_paths.R"
wd="/hpc/hub_oudenaarden/jyeung/code_for_analysis/scchic"
# outdir="/tmp"  # can change it
outdir=$1
[[ ! -d $1 ]] && echo "$1 not found, exiting" && exit 1
cd $wd

Rscript $rs $outdir
