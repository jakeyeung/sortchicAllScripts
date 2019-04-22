#!/bin/sh
# Jake Yeung
# 5-run_without_seed.sh
#  
# 2019-04-11

mark="H3K9me3"
Erds="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permutation_analysis_H3K9me3/mara_input/exprs_mat.rds"
Nrds="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permutation_analysis_H3K9me3/mara_input/sitecounts_mat.rds"
lambda="FALSE"
# indx="none"
# by="none"

indx=0
by="row"
marascript="/home/hub_oudenaarden/jyeung/projects/from_PhD/ridge-regression/run_ridge_regression2_permute.R"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permutation_analysis_${mark}/mara_output"

[[ ! -e $marascript ]] && echo "$marascript not found, exiting" && exit 1
[[ ! -e $Erds ]] && echo "$Erds not found, exiting" && exit 1
[[ ! -e $Nrds ]] && echo "$Nrds not found, exiting" && exit 1

outdir=$outmain/seed_${by}_${indx}
[[ ! -d $outdir ]] && mkdir $outdir
cd $outdir; Rscript $marascript $Erds $Nrds $lambda $indx $by

