#!/bin/sh
# Jake Yeung
# run.2-CCA_full_data_repressive_k9me3.sh
#  
# 2022-07-21

keeptop=500
jmem='32G'
jtime='1:00:00'
# rs="/nfs/scistore12/hpcgrp/jyeung/projects/scchic/scripts/revision_scripts/revisions2/2-CCA_full_data_repressive_k9me3_cluster_specific_bins_keeptop.R"
rs="/nfs/scistore12/hpcgrp/jyeung/projects/scchic/scripts/revision_scripts/revisions2/2-CCA_full_data_repressive_k9me3_cluster_specific_bins_keeptop_bymark_with_TSS.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

jmarkref="k9me3"
jmarks="k4me1 k4me3 k27me3"


jfactors="-1 1"

for jfactor in $jfactors; do
		outdir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/with_k9me3_cluster_specific_bins_keeptop_${keeptop}_bymark_factor_${jfactor}_with_TSS"
		[[ ! -d $outdir ]] && mkdir $outdir
		for jmark in $jmarks; do
		  echo $jmark

		  BNAME="${outdir}/sbatch_output_cca_${jmark}_with_TSS"
		  DBASE=$(dirname "${BNAME}")
		  [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
		  
		  module load R/4.1.2
		  cmd="Rscript $rs -refmark $jmarkref -mark $jmark -outdir $outdir -repressfactor $jfactor"
		  sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=${jmark}_${jfactor} --wrap "$cmd"
		done
done
