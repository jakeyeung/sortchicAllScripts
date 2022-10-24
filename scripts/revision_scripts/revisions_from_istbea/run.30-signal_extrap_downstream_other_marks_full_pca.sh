#!/bin/sh
# Jake Yeung
# run.30-signal_extrap_downstream_other_marks_full_pca.sh
#  
# 2022-05-05

jmem='32G'
jtime='12:00:00'

jmarks="k4me1 k4me3 k27me3 k9me3"

rs="/nfs/scistore12/hpcgrp/jyeung/projects/scchic/scripts/revision_scripts/revisions_from_istbea/30-signal_extrap_downstream_other_marks_full_pca.R"
# outdir="/nfs/scistore12/hpcgrp/jyeung/projects/scchic/scripts/revision_scripts/revisions_from_istbea/30-signal_extrap_downstream_other_marks_full_pca.R"
outdir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_downstream/with_umap_arrows"
[[ ! -d $outdir ]] && mkdir $outdir
ncores=4
for jmark in $jmarks; do
  BNAME=${outdir}/downstream_${jmark}_fullpca.log
  DBASE=$(dirname "${BNAME}")
  [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
  module load R/4.1.2
  cmd="Rscript $rs -mark $jmark -outdir $outdir"
  sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=${ncores} --job-name=${jmark} --wrap "$cmd"
done
