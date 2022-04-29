#!/bin/sh
# Jake Yeung
# 4-run_LDA_doubles.sh
#  
# 2021-11-21

jmem='16G'
jtime='48:00:00'

rs="/nfs/scistore12/hpcgrp/jyeung/projects/scChIX/utils/run_LDA_model2.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

ncores=1
topics="30"
topicsName=`echo $topics | sed 's/,/_/g'`
# binarize="FALSE"

outmain0="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k4me3_cleaned"
# outmain0"/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mats_k4me3_remove_bad_louvains_dynbins_TSS"
[[ ! -d $outmain0 ]] && mkdir $outmain0

inmain="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mats_k4me3_remove_bad_louvains_dynbins_TSS"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1

for inf in `ls -d ${inmain}/*.rds`; do
    prefix=$(basename $inf)
    prefix=${prefix%.*}
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    echo $inf
    outmain="${outmain0}"
    [[ ! -d $outmain ]] && mkdir $outmain
    [[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

    bname="${prefix}"

    outdir="${outmain}/lda_outputs.${bname}"
    [[ -d $outdir ]] && echo "$outdir found, continuing" && continue
    [[ ! -d $outdir ]] && mkdir -p $outdir

    BNAME=$outdir/$bname
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

	module load R/4.1.2
    cmd="module add R/4.1.2; Rscript $rs $inf $outdir --topics $topics --projname $bname"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=LDA_${bname} --wrap "$cmd"
done
