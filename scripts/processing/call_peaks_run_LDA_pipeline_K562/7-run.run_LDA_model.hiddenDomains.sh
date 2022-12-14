#!/bin/sh
# Jake Yeung
# 2b-run.run_LDA_model.H3K27me3.BroadPeaks.sh
#  
# 2018-12-29

jmem='96G'
jtime='24:00:00'

workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"

cd $workdir

rs="scripts/processing/lib/run_LDA_model.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

marks="H3K27me3 H3K4me3"
mindist="1000"
cell="K562"

K=20  # kind of useless parameter
topics="10,15,20,25"
topicsName=`echo $topics | sed 's/,/_/g'`
tunemodels="TRUE"
binarize="FALSE"
ncores=4
meth="hiddenDomains"

for mark in $marks; do
    inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/count_mat_${cell}/count_mats.${meth}.${mindist}/PZ-${cell}-${mark}.merged.NoCountThres.Robj"
    # meanmax="0.15"  # about 0.32 for pvalcutoff 0.5, 0.15 for pvalcutoff 0.3
    meanmax="1"  # meanmax = 1 keeps the Erdr1 and Mid1 genes, which are probably real? But there may be some weird peaks as well that are skewing results.
    # cellmin="1000"
    # cellmax="50000"
    # relax assumptions to capture more H3K4me3 and H3K9me3 cells?
    cellmin="100"
    cellmax="500000"
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/LDA_out_${cell}/ldaAnalysis${meth}_${mindist}/lda_outputs.${meth}.meanfilt_${meanmax}.cellmin_${cellmin}.cellmax_${cellmax}.binarize.${binarize}"
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
