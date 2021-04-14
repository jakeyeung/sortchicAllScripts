#!/bin/sh
# Jake Yeung
# 3b-run_unmixing_downstream_split_reads.sh
#  
# 2020-02-06

echo "Sleep 2 hours"
sleep 7200

# WRAP UP
while [[ `squeue -u jyeung | grep RunFits | wc -l` > 0 ]]; do
        echo "sleep for 60 seconds"
        sleep 60
done

# # WRAP UP
# while [[ `qstat | wc -l` > 0 ]]; do
#         echo "sleep for 60 seconds"
#         sleep 60
# done

jmem='32G'
jtime='3:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/memux_scripts/unmixing_downstream_split_reads_script.R"

mark1="H3K4me1"
mark2="H3K9me3"
markdbl="H3K4me1xH3K9me3"  # match lda but not input in Rscript

experi="mouse_spikein_BMround2all.dbl"
suffix="common_rows_match_dbl"

dblmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis"
# ldadir="/hpc/hub_oudenaarden/jyeung/data/dblchic/from_cluster/LDA_outputs/ldaAnalysisBins_BM_${experi}"
ldadir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all.dbl_common_rows"
inputdir="${dblmain}/double_staining_input/${experi}_${suffix}"
outputdir="${dblmain}/double_staining_output/${experi}_${suffix}"

# LDA: 
# /hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all.dbl_common_rows/lda_outputs.count_mat.H3K9me3.match_dbl.K-30.binarize.FALSE/ldaOut.count_mat.H3K9me3.match_dbl.K-30.Robj
fnamesuffix="match_dbl.K-30.Robj"
dnamesuffix="match_dbl.K-30.binarize.FALSE"

[[ ! -d $ldadir ]] && echo "$ldadir not found, exiting" && exit 1
[[ ! -d $inputdir ]] && echo "$inputdir not found, exiting" && exit 1
[[ ! -d $outputdir ]] && echo "$outputdir not found, exiting" && exit 1

    unmixingname="${experi}_${suffix}_clstr_by_louvain_${markdbl}.removeNA_FALSE"
    infdblinput="${inputdir}/${unmixingname}.RData"
    infdbloutput="${outputdir}/unmix_${unmixingname}.RData"

    [[ ! -e $infdblinput ]] && echo "$infdblinput not found, exiting" && exit 1
    [[ ! -e $infdbloutput ]] && echo "$infdbloutput not found, exiting" && exit 1

    inflda1="${ldadir}/lda_outputs.count_mat.${mark1}.${dnamesuffix}/ldaOut.count_mat.${mark1}.${fnamesuffix}"
    inflda2="${ldadir}/lda_outputs.count_mat.${mark2}.${dnamesuffix}/ldaOut.count_mat.${mark2}.${fnamesuffix}"
    infldadbl="${ldadir}/lda_outputs.count_mat.${markdbl}.${dnamesuffix}/ldaOut.count_mat.${markdbl}.${fnamesuffix}"

    [[ ! -e $inflda1 ]] && echo "$inflda1 not found, exiting" && exit 1
    [[ ! -e $inflda2 ]] && echo "$inflda2 not found, exiting" && exit 1
    [[ ! -e $infldadbl ]] && echo "$infldadbl not found, exiting" && exit 1

    outdir="${dblmain}/double_staining_output_downstream/${experi}_${suffix}/MF_BM_${experi}_${suffix}_${mark1}-${mark2}_split_reads"
    [[ ! -d $outdir ]] && mkdir -p $outdir

    bname="MF_${experi}_${suffix}_clstr_by_louvain.${mark1}x${mark2}"

    outprefix="${outdir}/${bname}"

    BNAME=$outprefix.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    # inputs are keywords
    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -mark1 $mark1 -mark2 $mark2 -inf_dbl_input $infdblinput -inf_dbl_output $infdbloutput -inf_mark1_lda $inflda1 -inf_mark2_lda $inflda2 -inf_dblmark_lda $infldadbl -outprefix $outprefix" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N unmix_ds_${markdbl}
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs -mark1 $mark1 -mark2 $mark2 -inf_dbl_input $infdblinput -inf_dbl_output $infdbloutput -inf_mark1_lda $inflda1 -inf_mark2_lda $inflda2 -inf_dblmark_lda $infldadbl -outprefix $outprefix"
    # | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N unmix_ds_${markdbl}
    echo $cmd
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=unmix_ds_${markdbl} --wrap "$cmd"


