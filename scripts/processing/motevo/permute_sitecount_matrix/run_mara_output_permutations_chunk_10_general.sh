#!/bin/sh
# Jake Yeung
# 2-run_mara_output_permutations.sh
# Run a chunk of 10
# 2019-04-01

deci=$1
[[ $deci != [0-9]* ]] && echo "Must be integer: $deci" && exit 1
by=$2
mark=$3
Erds=$4
Nrds=$5
outmain=$6

indxstart=$(expr $deci \* 10 - 9)
indxend=$(expr $deci \* 10 - 9 + 9) 

echo "$indxstart - $indxend"

marascript="/home/hub_oudenaarden/jyeung/projects/from_PhD/ridge-regression/run_ridge_regression2_permute.R"
[[ ! -e $marascript ]] && echo "$marascript not found, exiting" && exit 1

# run 10 jobs
# Erds="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permutation_analysis_${mark}/mara_input/exprs_mat.rds"
# Nrds="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permutation_analysis_${mark}/mara_input/sitecounts_mat.rds"

[[ ! -e $Erds ]] && echo "$Erds not found, exiting" && exit 1
[[ ! -e $Nrds ]] && echo "$Nrds not found, exiting" && exit 1

lambda="FALSE"
# outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permutation_analysis_${mark}/mara_output"

[[ ! -d $outmain ]] && mkdir -p $outmain

# by="column"
for indx in $(seq $indxstart $indxend); do
    outdir=$outmain/seed_${by}_${indx}
    [[ ! -d $outdir ]] && mkdir $outdir
    cd $outdir; Rscript $marascript $Erds $Nrds $lambda $indx $by&
    # cd $outdir; Rscript $marascript $Erds $Nrds $lambda $indx $by
    # echo "Debug. Exiting"
    # exit 0
done
wait
