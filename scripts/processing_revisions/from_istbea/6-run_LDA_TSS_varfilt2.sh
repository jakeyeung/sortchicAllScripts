#!/bin/sh
# Jake Yeung
# 4-run_LDA_doubles.sh
#  
# 2021-11-21

jmem='16G'
jtime='96:00:00'
# jtime='1:00:00'

# conda deactivate

rs="/nfs/scistore12/hpcgrp/jyeung/projects/scChIX/utils/run_LDA_model2.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

ncores=1
topics="30"
topicsName=`echo $topics | sed 's/,/_/g'`
# binarize="FALSE"

outmain0="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_bugfixed_TSS"
[[ ! -d $outmain0 ]] && mkdir $outmain0

# jnames="k4me1 k4me3 k27me3 k9me3"
# jnames="BM_k4me1"
# jname="BM_k4me1"

# inmain="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mat_TSS"
# inmain="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mat_TSS_rep2rep3only"
inmain="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mat_TSS_varfilt2"
[[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
outmain="${outmain0}/ldaAnalysis_varfilt_again"
[[ ! -d $outmain ]] && mkdir $outmain

for inf in `ls -d ${inmain}/*.rds`; do
    prefix=$(basename $inf)
    prefix=${prefix%.*}
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    echo $inf
    [[ ! -d $outmain ]] && mkdir $outmain
    [[ ! -d $outmain ]] && echo "$outmain not found, exiting" && exit 1
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    
    bname="${prefix}"
    
    # outdir="${outmain}/lda_outputs.${bname}"
	outdir=${outmain}
    # [[ -d $outdir ]] && echo "$outdir found, continuing" && continue
    # [[ ! -d $outdir ]] && mkdir -p $outdir
    
    BNAME=$outdir/$bname
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    
    module load R/4.1.2
    cmd="module load R/4.1.2; which R; Rscript $rs $inf $outdir --topics $topics --projname $bname"
    echo $cmd
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=LDA_${bname} --wrap "$cmd"
done

