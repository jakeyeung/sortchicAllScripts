#!/bin/sh
# Jake Yeung
# 1-project_new_samples_on_LDA.H3K27me3.StemCells.sh
#  
# 2019-10-24



jmem='24G'
jtime='24:00:00'

rs="/home/hub_oudenaarden/jyeung/projects/scchic-functions/scripts/processing_scripts/project_new_samples_on_LDA_bin.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# jmarkref="H3K27me3"
jmarkref="H3K4me1"
jmarks="K4m1 K27m3 K4m1_K27m3"
# markdbl="K4m1_K27m3"

# outdir="/hpc/hub_oudenaarden/jyeung/data/dblchic/from_cluster/projections_LDA_outputs/${experi}_${experilong}"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_projections"
[[ ! -d $outdir ]] && mkdir $outdir
[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1

ldadir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_BMAllMerged.2020-02-15.from_hiddendomains.NewCountFilters/lda_outputs.merged.${jmarkref}.minlength_1000.cutoff_analysis.merged.withchr.annotated.NewCountFilters.K-30.binarize.FALSE"
ldaname="ldaOut.merged.${jmarkref}.minlength_1000.cutoff_analysis.merged.withchr.annotated.NewCountFilters.K-30.Robj"
inlda="${ldadir}/${ldaname}"
[[ ! -e $inlda ]] && echo "$inlda not found, exiting" && exit 1

inmatdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_mat_B6_from_chix/countTables.unfixed"
for jmark in $jmarks; do
    inmat="${inmatdir}/all_BM_${jmark}_200119.mq_40.on_${jmarkref}_peaks.rds"
    outf="${outdir}/project_unmixed_${jmark}.project_on_AllMerged_${jmarkref}.RData"

    BNAME=${outdir}/${jmark}_on_${jmarkref}.qsub
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    echo ". /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inlda $inmat $outf" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 1 -m beas -M j.yeung@hubrecht.eu -N ${jmark}_project_LDA


done

