#!/bin/sh
# Jake Yeung
# 1-project_new_samples_on_LDA.H3K27me3.StemCells.sh
#  
# 2019-10-24

# echo "Sleep 12000"
# sleep 12000

# WRAP UP
while [[ `squeue -u jyeung | grep RunFits | wc -l` > 0 ]]; do
        echo "sleep for 180 seconds"
        sleep 180
done

echo "Sleep another 600 to be sure"
sleep 600

# WRAP UP
while [[ `squeue -u jyeung | grep unmix_ds | wc -l` > 0 ]]; do
        echo "sleep for 180 seconds"
        sleep 180
done

jmem='24G'
jtime='24:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/project_new_samples_on_LDA_bin.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# jmarks="K27m3 K9m3"
mark1="H3K4me1"
mark2="H3K9me3"
jmarks="H3K4me1 H3K9me3"

# prefix="EtOH_NoTcells_VarFilt_pass2_autosomesOnly"
# prefix="EtOH_NoTcells_VarFilt_pass2_autosomesOnly.maxcountsfilt"
# ldaAnalysisBins_mouse_spikein_BMround2all.dbl_common_rows

# prefix="mouse_spikein_BMround2all.dbl_common_rows"
# suffix="match_dbl"

prefix="SetupObjs_AllMerged"
# suffix="UnionRows"
suffix="UnionRows_KeepAllCells_FewerRepressedClusters_MergeActiveClusters3"
# experi=${prefix}_${suffix}

markdbl="H3K4me1xH3K9me3"
jsuffix="match_dbl"

# experilong="BM_EtOH_K27m3-K9m3_split_reads"

# outdir="/hpc/hub_oudenaarden/jyeung/data/dblchic/from_cluster/projections_LDA_outputs/afterUnmixing_${experi}_${markdbl}"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/projects_after_unmixing.${markdbl}"
[[ ! -d $outmain ]] && mkdir $outmain
outdir=${outmain}/${prefix}_${suffix}
[[ ! -d $outdir ]] && mkdir $outdir
[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1


# ldadir="/hpc/hub_oudenaarden/jyeung/data/dblchic/from_cluster/LDA_outputs/ldaAnalysisBins_BM_${prefix}"
# ldadir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_${prefix}"
ldadir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld"
# dnamesuffix="KeepTopBins_500.KeepAllPlates.K-30.binarize.FALSE"
# fnamesuffix="KeepTopBins_500.KeepAllPlates.K-30.Robj"

for jmark in $jmarks; do
    dname="lda_outputs.count_mat_old_merged_with_new.${jmark}.K-30.binarize.FALSE"
    # dname="lda_outputs.count_mat.${jmark}.${suffix}.K-30.binarize.FALSE"
    fname="ldaOut.count_mat_old_merged_with_new.${jmark}.K-30.Robj"
    # fname="ldaOut.count_mat.${jmark}.${suffix}.K-30.Robj"
    inlda=${ldadir}/${dname}/${fname}

    # inmat="/hpc/hub_oudenaarden/jyeung/data/dblchic/double_staining_output_downstream/${prefix}_${suffix}/MF_BM_${prefix}_${suffix}_K27m3-K9m3_split_reads/MF_${prefix}_${suffix}_clstr_by_louvain.K27m3-K9m3-unmixed_mat.${jmark}.rds"
    inmat="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/H3K4me1_H3K9me3_analysis/double_staining_output_downstream/${prefix}_${suffix}/MF_BM_${prefix}_${suffix}_${mark1}-${mark2}_split_reads/MF_${prefix}_${suffix}_clstr_by_louvain.${markdbl}-unmixed_mat.${jmark}.rds"

    [[ ! -e $inmat ]] && echo "$inmat not found, exiting" && exit 1
    [[ ! -e $inlda ]] && echo "$inlda not found, exiting" && exit 1

    outf="${outdir}/project_unmixed_${jmark}.RData"

    BNAME=${outdir}/${jmark}_qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    # echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inlda $inmat $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N ${jmark}_project_LDA
    cmd=". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inlda $inmat $outf"
    sbatch --time=$jtime --mem-per-cpu=$jmem --output=${BNAME}_%j.log --ntasks=1 --nodes=1 --job-name=project_${jmark} --wrap "$cmd"

done

