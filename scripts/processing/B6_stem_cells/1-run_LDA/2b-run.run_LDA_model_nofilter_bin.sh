#!/bin/sh
# Jake Yeung
# 2b-run.run_LDA_model.H3K27me3.BroadPeaks.sh
#  
# 2018-12-29

# if not binarize
# jmem='64G'
# jtime='24:00:00'
# if binarize


jmem='48G'
jtime='51:00:00'

workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"

cd $workdir

rs="scripts/processing/lib/run_LDA_model.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

marks="H3K4me1"

K=20  # kind of useless parameter
ncores=3
topics="30,35,50"
topicsName=`echo $topics | sed 's/,/_/g'`
tunemodels="TRUE"
binarize="TRUE"
cellmin="0"
cellmax="999999"
# relax assumptions to capture more H3K4me3 and H3K9me3 cells?
meanmax="999"  # meanmax = 1 keeps the Erdr1 and Mid1 genes, which are probably real? But there may be some weird peaks as well that are skewing results.

prefix="PZ-Bl6-BM-StemCells"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA_${prefix}"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

mindist="$prefix"

for mark in $marks; do
    for inf in `ls -d $indir/${prefix}_${mark}*.RData`; do
        outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_${mindist}/lda_outputs.meanfilt_${meanmax}.cellmin_${cellmin}.cellmax_${cellmax}.binarize.${binarize}.no_filt"
        [[ ! -d $outdir ]] && mkdir -p $outdir
        [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
        [[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1
        bname=$(basename $inf)
        bname=${bname%%.*}.CountThres0.K-${topicsName}
        BNAME=$outdir/$bname
        DBASE=$(dirname "${BNAME}")
        [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
        echo "cd $workdir; /hpc/hub_oudenaarden/jyeung/software/anaconda3/envs/R3.6/bin/Rscript $rs $inf $outdir $K $topics $tunemodels $meanmax $cellmin $cellmax $binarize $bname" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu
    done
done
