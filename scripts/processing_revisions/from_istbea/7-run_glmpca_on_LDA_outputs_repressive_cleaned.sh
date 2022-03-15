#!/bin/sh
# Jake Yeung
# 1-run_glmpca_using_LDA_input.sh
#  
# 2020-11-17

jmem='32G'
jtime='24:00:00'

rs="/nfs/scistore12/hpcgrp/jyeung/git_repos/scchic-functions/scripts/processing_scripts/run_GLMPCA_with_LDA_init_spikeins_plate.R"

jmarks="k9me3"
platename="batch"
szname="none"

bincutoff=0
binskeepvec="0"
nitervec="1000"

outdir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/glmpca_outputs/repressive_cleaned"
[[ ! -d $outdir ]] && mkdir $outdir
annotdir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata/repressive_cleaned"
ldadir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_repressive_cleaned_from_jupyter/ldaAnalysis_fripfilt_varfilt_binfilt"

for binskeep in $binskeepvec; do
    for niter in $nitervec; do
        echo $binskeep
        echo $niter

        for jmark in $jmarks; do
			# jbase="count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.${jmark}.2022-02-13"
			jbase="count_mat_cleaned_dynbins.${jmark}.2022-02-16"
			infpeaks="${ldadir}/lda_outputs.${jbase}/ldaOut.${jbase}.Robj"
			[[ ! -e $infpeaks ]] && echo "$infpeaks not found, exiting" && exit 1

            outbase=${outdir}/glmpca.${jmark}.bincutoff_${bincutoff}.binskeep_${binskeep}.platename_${platename}.szname_${szname}.niter_${niter}
            # check output doesnt exist
            outcheck="${outbase}.RData"
            [[ -e $outcheck ]] && echo "$outcheck found, continuing" && continue
            BNAME=${outdir}/glmpca.${jmark}.bincutoff_${bincutoff}.binskeep_${binskeep}.platename_${platename}.szname_${szname}.niter_${niter}
            DBASE=$(dirname "${BNAME}")
            [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

			# annotname="metadata_plate_experi_batch.${jmark}.2022-01-27.txt"
			annotname="metadata_plate_experi_batch.cleaned.dynamicbins.${jmark}.2022-02-17.txt"
            annotf=${annotdir}/${annotname}
			[[ ! -e $annotf ]] && echo "$annotf not found, exiting" && exit 1

			module load R/4.1.2
            cmd="Rscript $rs -infpeaks $infpeaks -infmeta $annotf -outbase $outbase -platecname $platename -sizefactorcname $szname -bincutoff $bincutoff -binskeep $binskeep -niter $niter"
            sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${jmark} --wrap "$cmd"
        done
    done
done


