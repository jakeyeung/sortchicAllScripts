#!/bin/sh
# Jake Yeung
# run.23-correct_batch_effect_k27me3_TSS.sh
#  
# 2022-04-24

jmem='72G'
jtime='24:00:00'

BNAME="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/correct_batch_k27me3_TSS.sbatch_log"
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

rs="/nfs/scistore12/hpcgrp/jyeung/projects/scchic/scripts/revision_scripts/revisions_from_istbea/23-correct_batch_effect_k27me3_TSS.R"

module load R/4.1.2
cmd="Rscript $rs"
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=4 --job-name=BatchK27me3 --wrap "$cmd"

