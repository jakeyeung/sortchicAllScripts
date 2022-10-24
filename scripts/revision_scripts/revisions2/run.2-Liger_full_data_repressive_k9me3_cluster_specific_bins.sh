#!/bin/sh
# Jake Yeung
# run.2-CCA_full_data_repressive_k9me3.sh
#  
# 2022-07-21

jmem='32G'
jtime='3:00:00'
# rs="/nfs/scistore12/hpcgrp/jyeung/projects/scchic/scripts/revision_scripts/revisions2/2-CCA_full_data_repressive_k9me3_cluster_specific_bins.R"
rs="/nfs/scistore12/hpcgrp/jyeung/projects/scchic/scripts/revision_scripts/revisions2/2-Liger_full_data_repressive_k9me3_cluster_specific_bins_keeptop_bymark.R"

# BNAME="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/with_k9me3/sbatch_output_cca"
outmain="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/liger_outputs/k9me3"
# BNAME="${outmain}/sbatch_output_liger"
# DBASE=$(dirname "${BNAME}")
# [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

jmarkref="k9me3"
jmarks="k4me1 k4me3 k27me3"
jfactors="-1 1"

keeptop=500
for jfactor in $jfactors; do
  outdir="${outmain}/with_k9me3_cluster_specific_bins_keeptop_${keeptop}_bymark_factor_${jfactor}"
  [[ ! -d $outdir ]] && mkdir $outdir
  for jmark in $jmarks; do
    echo $jmark
  
    BNAME="${outdir}/sbatch_output_liger_${jmark}"
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
  
    module load R/4.1.2
    module load hdf5/1.13.0
    cmd="Rscript $rs -refmark $jmarkref -mark $jmark -outdir $outdir -repressfactor $jfactor"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=${jmark}_${jfactor} --wrap "$cmd"
  done
done

# module load R/4.1.2
# module load hdf5/1.13.0
# cmd="Rscript $rs"
# sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=cca --wrap "$cmd"
