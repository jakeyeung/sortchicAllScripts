#!/bin/sh
# Jake Yeung
# run.make_fig2_server_paths.sh
#  
# 2019-05-25

outdir=$1
[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1

rs="/hpc/hub_oudenaarden/jyeung/code_for_analysis/scchic/scripts/scripts_analysis/B6_analysis/make_figs/make_fig2_server_paths.R"
wd="/hpc/hub_oudenaarden/jyeung/code_for_analysis/scchic"

cd $wd
Rscript $rs $outdir
