#!/bin/sh
# Jake Yeung
# run.zinbwave_and_deseq_on_pseudotime.sh
#  
# 2019-01-02

jmem='96G'
jtime='48:00:00'
BNAME="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_tmp/zinbwave_run"
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/scripts_analysis/zinb_deseq/zinbwave_and_deseq_on_pseudotime.R"

echo "Rscript $rs" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu
