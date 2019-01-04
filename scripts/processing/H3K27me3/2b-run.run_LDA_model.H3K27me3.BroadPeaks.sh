#!/bin/sh
# Jake Yeung
# 2b-run.run_LDA_model.H3K27me3.BroadPeaks.sh
#  
# 2018-12-29

jmem='64G'
jtime='24:00:00'

workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"

cd $workdir

rs="scripts/processing/run_LDA_model.R"

inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats/PZ-BM-H3K27me3.merged.NoCountThres.Robj"

dist="1000_NoM_binarize"
meanmax="1"
cellmin="1000"
cellmax="50000"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/lda_outputs.meanfilt_${meanmax}.merge_${dist}.cellmin_${cellmin}.cellmax_${cellmax}"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/lda_analysis.h3k27me3.broadpeaks/lda_outputs.meanfilt_${meanmax}.merge_${dist}.cellmin_${cellmin}.cellmax_${cellmax}"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/lda_analysis.h3k27me3.dropbox/lda_output"
[[ ! -d $outdir ]] && mkdir $outdir


# outf=$outdir/$bname.lda_out.RData
# outftune=$outdir/$bname.lda_tune.RData


[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1

K=20
topics="10,15,20,25,30,40,50"
tunemodels="TRUE"
binarize="TRUE"
ncores=7

bname=$(basename $inf)
bname=${bname%%.*}.CountThres0.K-${K}
BNAME=$outdir/$bname
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

echo "cd $workdir; Rscript $rs $inf $outdir $K $topics $tunemodels $meanmax $cellmin $cellmax $binarize" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded $ncores -m beas -M j.yeung@hubrecht.eu
