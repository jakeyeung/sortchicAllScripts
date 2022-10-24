#!/bin/sh
# Jake Yeung
# run.2-CCA_full_data_repressive_k9me3.sh
#  
# 2022-07-21

jmem='64G'
jtime='6:00:00'
# rs="/nfs/scistore12/hpcgrp/jyeung/projects/scchic/scripts/revision_scripts/revisions2/2-CCA_full_data_repressive_k9me3_cluster_specific_bins.R"
rs="/nfs/scistore12/hpcgrp/jyeung/projects/scchic/scripts/revision_scripts/revisions2/2-Liger_full_data_repressive.R"

# BNAME="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/with_k9me3/sbatch_output_cca"
BNAME="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/liger_outputs/sbatch_liger"
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1


# cmd="module load R/4.1.2; Rscript $rs"
module load R/4.1.2
module load hdf5/1.13.0
cmd="Rscript $rs"
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=liger --wrap "$cmd"
