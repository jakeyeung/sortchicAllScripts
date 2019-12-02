#!/bin/sh
# Jake Yeung
# 1-run.run_LDA_model_ZF.sh
# 2019-11-03

# if not binarize
# jmem='64G'
# jtime='24:00:00'
# if binarize

jmem='48G'
jtime='24:00:00'

workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"

cd $workdir

rs="scripts/processing/lib/run_LDA_model.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

# marks="H3K36me3"
# marks="H3K27me3"
# marks="H3K4me3"

K=20  # kind of useless parameter
ncores=3
topics="30,35,50"
topicsName=`echo $topics | sed 's/,/_/g'`
tunemodels="TRUE"
binarize="FALSE"
cellmin="0"
cellmax="999999"
# relax assumptions to capture more H3K4me3 and H3K9me3 cells?
meanmax="999"  # meanmax = 1 keeps the Erdr1 and Mid1 genes, which are probably real? But there may be some weird peaks as well that are skewing results.

prefix="ZFbonemarrow"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA_${prefix}"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1

outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_${prefix}"

jdate="2019-11-04"

for inf in `ls -d $indir/ZF-*.${jdate}.RData`; do
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    bname=$(basename $inf)
    bname=${bname%%.*}.CountThres0.K-${topicsName}

    outdir="${outmain}/lda_outputs.${bname}.binarize.${binarize}.${jdate}"
    [[ ! -d $outdir ]] && mkdir -p $outdir

    BNAME=$outdir/$bname
    DBASE=$(dirname "${BNAME}")
    [[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

    echo "cd $workdir; /hpc/hub_oudenaarden/jyeung/software/anaconda3/envs/R3.6/bin/Rscript $rs $inf $outdir $K $topics $tunemodels $meanmax $cellmin $cellmax $binarize $bname" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu -N $bname
done
