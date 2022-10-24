#!/bin/sh
# Jake Yeung
# run.5-summarize_signal_different_bin_sets.sh
#  
# 2022-07-25

jmem='64G'
jtime='12:00:00'

BNAME="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/counts_different_bins/summarize_signal"
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

module load R/4.1.2
rs="/nfs/scistore12/hpcgrp/jyeung/projects/scchic/scripts/revision_scripts/revisions2/5-summarize_signal_different_bin_sets.R"
cmd="Rscript $rs"
sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=Sum --wrap "$cmd"

