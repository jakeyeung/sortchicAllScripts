#!/bin/sh
# Jake Yeung
# run.2-CCA_full_data_repressive_k9me3.sh
#  
# 2022-07-21

jmem='96G'
jtime='3:00:00'
rs="/nfs/scistore12/hpcgrp/jyeung/projects/scchic/scripts/revision_scripts/revisions2/2-CCA_full_data_repressive_k9me3.R"

BNAME="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/with_k9me3/sbatch_output_cca"
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1


# cmd="module load R/4.1.2; Rscript $rs"
module load R/4.1.2
cmd="Rscript $rs"
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=cca --wrap "$cmd"
