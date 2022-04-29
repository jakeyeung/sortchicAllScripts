#!/bin/sh
# Jake Yeung
# 1-run_glmpca_using_LDA_input.sh
#  
# 2020-11-17

jmem='32G'
jtime='24:00:00'

# rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/run_GLMPCA_with_LDA_init_spikeins_plate.from_imputed.R"
rs="/nfs/scistore12/hpcgrp/jyeung/projects/scchic-functions/scripts/processing_scripts/run_GLMPCA_with_LDA_init_spikeins_plate.from_imputed.R"

jmarks="k27me3"
platename="batch"
szname="none"

nitervec="11"

# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/glmpca_outputs/BM_from_matadj"
outdir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/glmpca_outputs/batch_corrected_k27me3"
[[ ! -d $outdir ]] && mkdir $outdir

annotdir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths2"

indirimputed="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections"
indircounts="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections"

# jsuffix="_old_to_new"
# jsuffix2="TES"

for niter in $nitervec; do
    echo $niter
    for jmark in $jmarks; do
        fnameimputed="mat_wide_${jmark}_batch_corrected.2022-04-19.rds"
        fnamecounts="mat_wide_${jmark}_raw_counts.2022-04-19.rds"
        infimputed="${indirimputed}/$fnameimputed"
        infcounts="${indircounts}/$fnamecounts"

        outbase=${outdir}/glmpca.${jmark}.from_matadj.platename_${platename}.szname_${szname}.niter_${niter}
        outcheck="${outbase}.RData"
        [[ -e $outcheck ]] && echo "$outcheck found, continuing" && continue
        BNAME=${outdir}/glmpca.${jmark}.from_matadj.platename_${platename}.szname_${szname}.niter_${niter}.sbatch
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

        # annotname="count_mat_from_sitecount_mat.${jmark}.filtNAcells_allbins.from_same_annot_file.metadata.2020-12-27.txt"
		annotname="metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.${jmark}.txt"
        annotf=${annotdir}/${annotname}
		module load R/4.1.2
        cmd="Rscript $rs -infimputed $infimputed -infcounts $infcounts -infmeta $annotf -outbase $outbase -platecname $platename -sizefactorcname $szname -niter $niter -jntopics 30"
        sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${jmark} --wrap "$cmd"
    done
done
