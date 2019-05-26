#!/bin/sh
# Jake Yeung
# run.make_fig3_server_paths.sh
# Run it 
# 2019-05-26

bs="/hpc/hub_oudenaarden/jyeung/code_for_analysis/scchic/scripts/scripts_analysis/B6_analysis/make_figs/make_fig3_server_paths.R"
[[ ! -e $bs ]] && echo "$bs not found, exiting" && exit 1

outdir="/hpc/hub_oudenaarden/jyeung/tmp"

[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1

wd="/hpc/hub_oudenaarden/jyeung/code_for_analysis/scchic"

cd $wd

Rscript $bs $outdir
