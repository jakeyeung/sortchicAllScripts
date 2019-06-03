#!/bin/sh
# Jake Yeung
# run.make_fig4a_umap_with_umi_counts.sh
#  
# 2019-05-31

rs="/hpc/hub_oudenaarden/jyeung/code_for_analysis/scchic/scripts/scripts_analysis/B6_analysis/make_figs/make_fig4a_umap_with_umi_counts.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
wd="/hpc/hub_oudenaarden/jyeung/code_for_analysis/scchic"
[[ ! -d $wd ]] && echo "$wd not found, exiting" && exit 1
outdir=$1
[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1
cd $wd
Rscript $rs $outdir
