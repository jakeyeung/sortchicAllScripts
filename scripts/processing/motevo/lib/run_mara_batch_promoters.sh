#!/bin/sh
# Jake Yeung
# run_mara_batch_promoters.sh
# Run MARA given E, output, N, lambda is optional
# 2018-11-29

# wd="/home/hub_oudenaarden/jyeung/projects/from_PhD/ridge-regression"
exprs_mat=$1
out=$2
sitecountmat=$3
lambda=$4  # can set to nothing to automatically find lambda through cross validation

[[ ! -e $sitecountmat ]] && echo "$sitecountmat not found, exiting" && exit 1
[[ ! -e $exprs_mat ]] && echo "$exprs_mat not found" && exit 1


dirname=$(basename $exprs_mat)
dirname=${dirname%.*}
outdir=$out/"$dirname"

marascript="/home/hub_oudenaarden/jyeung/projects/from_PhD/ridge-regression/run_ridge_regression2.R"

[[ ! -e $marascript ]] && echo "$marascript not found, exiting" && exit 1
[[ ! -d $outdir ]] && mkdir $outdir
[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1

N=$sitecountmat
E=$exprs_mat

cd $outdir

Rscript $marascript $E $N $lambda
