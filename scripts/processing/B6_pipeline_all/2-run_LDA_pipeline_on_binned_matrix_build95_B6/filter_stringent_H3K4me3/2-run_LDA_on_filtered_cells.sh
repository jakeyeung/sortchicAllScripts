#!/bin/sh
# Jake Yeung
# 2-run_LDA_on_filtered_cells.sh
#  
# 2019-06-03

jmem='48G'
jtime='51:00:00'

workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"

cd $workdir

rs="scripts/processing/lib/run_LDA_model.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

marks="H3K4me3"
mindist="BinFiltCellFilt2"
# cell="BM"

K=20  # kind of useless parameter
ncores=4
topics="25,30,35,50"
topicsName=`echo $topics | sed 's/,/_/g'`
tunemodels="TRUE"
binarize="TRUE"
cellmin="0"
cellmax="999999"
# relax assumptions to capture more H3K4me3 and H3K9me3 cells?
meanmax="999"  # meanmax = 1 keeps the Erdr1 and Mid1 genes, which are probably real? But there may be some weird peaks as well that are skewing results.

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA/stringent_filter"
prefix="B6"

# marks="H3K27me3"
for mark in $marks; do
    # inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats_binned/${cell}-${mark}.AvO_filt.Robj"
    # inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats_binned/${cell}-${mark}.no_filt.Robj"
    inf=$indir/"${prefix}_${mark}_pcutoff_0.95_binfilt_cellfilt.2019-06-03.stringent_filter.RData"
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_${mindist}/lda_outputs.meanfilt_${meanmax}.cellmin_${cellmin}.cellmax_${cellmax}.binarize.${binarize}.no_filt.stringent_filter"
    [[ ! -d $outdir ]] && mkdir -p $outdir
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    [[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1

    bname=$(basename $inf)
    bname=${bname%%.*}.CountThres0.K-${topicsName}
    BNAME=$outdir/$bname
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
    # echo $bname
    echo "cd $workdir; Rscript $rs $inf $outdir $K $topics $tunemodels $meanmax $cellmin $cellmax $binarize $bname" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu
    # echo "cd $workdir; Rscript $rs $inf $outdir $K $topics $tunemodels $meanmax $cellmin $cellmax $binarize $bname"
done
