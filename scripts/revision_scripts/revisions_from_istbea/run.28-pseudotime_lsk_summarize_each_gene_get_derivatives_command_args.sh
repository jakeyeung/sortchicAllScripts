#!/bin/sh
# Jake Yeung
# run.28-pseudotime_cubic_spline_each_gene_get_derivatives.sh
#  
# 2022-05-04

jmem='64G'
jtime='12:00:00'
# rs="/nfs/scistore12/hpcgrp/jyeung/projects/scchic/scripts/revision_scripts/revisions_from_istbea/28-pseudotime_lsk_cubic_spline_each_gene_get_derivatives_command_args.R"
rs="/nfs/scistore12/hpcgrp/jyeung/projects/scchic/scripts/revision_scripts/revisions_from_istbea/28-pseudotime_lsk_summarize_each_gene_get_derivatives_command_args.R"

jmarks="k4me1 k4me3 k27me3 k9me3"
outdir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/gam_fits_lsk"
[[ ! -d $outdir ]] && mkdir $outdir

for jmark in $jmarks; do

	outrds=${outdir}/"gam_fits_lsk_dynamic_bins.${jmark}.fromargs.full_list.summarize_mean_mad_sd.2022-05-07.rds"

    BNAME=${outdir}/sbatch_${jmark}.out
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1


	module load R/4.1.2
	cmd="Rscript $rs -mark $jmark -outfile $outrds"
	sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --cpus-per-task=1 --job-name=${jmark} --wrap "$cmd"
done
