#!/bin/sh
# Jake Yeung
# 1-run_glmpca_using_LDA_input.sh
#  
# 2020-11-17

jmem='32G'
jtime='24:00:00'

rs="/nfs/scistore12/hpcgrp/jyeung/git_repos/scchic-functions/scripts/processing_scripts/run_GLMPCA_with_LDA_init_spikeins_plate.R"

jmarks="k4me1 k4me3 k27me3 k9me3"
# jmarks="k27me3"
platename="batch"
szname="none"

bincutoff=0
binskeepvec="0"
nitervec="500 1000"

outdir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/glmpca_outputs/varfilt"
[[ ! -d $outdir ]] && mkdir $outdir
# annotdir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata"
annotdir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata/varfilt"
ldadir="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt"

for binskeep in $binskeepvec; do
    for niter in $nitervec; do
        echo $binskeep
        echo $niter

        for jmark in $jmarks; do
			jbase="count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.${jmark}.2022-01-28"
			infpeaks="${ldadir}/lda_outputs.${jbase}/ldaOut.${jbase}.Robj"
			[[ ! -e $infpeaks ]] && echo "$infpeaks not found, exiting" && exit 1
			# infpeaks="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_bugfixed_dynamic_bins/ldaAnalysis_fripfilt_BM_${jmark}/lda_outputs.count_mat_merged_with_old.${jmark}.2022-01-26/ldaOut.count_mat_merged_with_old.${jmark}.2022-01-26.Robj"
			# infpeaks="/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.k9me3.2022-01-28/ldaOut.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.k9me3.2022-01-28.Robj"

            outbase=${outdir}/glmpca.${jmark}.bincutoff_${bincutoff}.binskeep_${binskeep}.platename_${platename}.szname_${szname}.niter_${niter}
            # check output doesnt exist
            outcheck="${outbase}.RData"
            [[ -e $outcheck ]] && echo "$outcheck found, continuing" && continue
            BNAME=${outdir}/glmpca.${jmark}.bincutoff_${bincutoff}.binskeep_${binskeep}.platename_${platename}.szname_${szname}.niter_${niter}
            DBASE=$(dirname "${BNAME}")
            [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

			# annotname="metadata_plate_experi_batch.${jmark}.2022-01-27.txt"
			annotname="metadata_plate_experi_batch.varfilt.dynamicbins.${jmark}.2022-01-29.txt"
            annotf=${annotdir}/${annotname}
			[[ ! -e $annotf ]] && echo "$annotf not found, exiting" && exit 1

			module load R/4.1.2
            cmd="Rscript $rs -infpeaks $infpeaks -infmeta $annotf -outbase $outbase -platecname $platename -sizefactorcname $szname -bincutoff $bincutoff -binskeep $binskeep -niter $niter"
            sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=${jmark} --wrap "$cmd"
        done
    done
done


