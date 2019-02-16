#!/bin/sh
# Jake Yeung
# 2b-run.run_LDA_model.H3K27me3.BroadPeaks.sh
#  
# 2018-12-29

# if not binarize
# jmem='96G'
# jtime='8:00:00'
# if binarize

jmem='64G'
jtime='48:00:00'

workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"

cd $workdir

rs="scripts/processing/lib/run_LDA_model.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

marks="H3K4me1 H3K27me3 H3K9me3 H3K4me3"
# marks="H3K27me3"
mindist="MetaCell"
cell="K562-G1"

K=20  # kind of useless parameter
ncores=4
topics="15,20,25,30"
topicsName=`echo $topics | sed 's/,/_/g'`
tunemodels="TRUE"
binarize="FALSE"
cellmin="100"
cellmax="500000"
# relax assumptions to capture more H3K4me3 and H3K9me3 cells?
meanmax="10"  # meanmax = 1 keeps the Erdr1 and Mid1 genes, which are probably real? But there may be some weird peaks as well that are skewing results.

# marks="H3K27me3"
for mark in $marks; do
    inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_mat_K562/count_mats_binned/${cell}-${mark}.no_filt.Robj"
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/LDA_out_K562/ldaAnalysisBins_${mindist}/lda_outputs.meanfilt_${meanmax}.cellmin_${cellmin}.cellmax_${cellmax}.binarize.${binarize}.no_filt"
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

done
