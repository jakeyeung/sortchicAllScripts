#!/bin/sh
# Jake Yeung
# 7-run.run_LDA_model.H3K27me3.sh
# Run model after Alexander regenerated binned counts for LDA
# Filter cells and peaks based on metacell output
# 2018-12-28

jmem='16G'
jtime='24:00:00'

workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"

cd $workdir

rs="scripts/processing/lib/run_LDA_model.R"

inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/lda_analysis.h3k27me3.dropbox/BM-H3K27me3.AvO_filt.Robj"

dist="Windows2"
meanmax="1000"
cellmin="NA"
cellmax="NA"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/lda_outputs.meanfilt_${meanmax}.merge_${dist}.cellmin_${cellmin}.cellmax_${cellmax}"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/lda_analysis.h3k27me3.dropbox/lda_outputs.meanfilt_${meanmax}.merge_${dist}.cellmin_${cellmin}.cellmax_${cellmax}"
# outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/lda_analysis.h3k27me3.dropbox/lda_output"
[[ ! -d $outdir ]] && mkdir $outdir


# outf=$outdir/$bname.lda_out.RData
# outftune=$outdir/$bname.lda_tune.RData


[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1

K=12
topics="10,15,20,50,75,100"
tunemodels="FALSE"

bname=$(basename $inf)
bname=${bname%%.*}.CountThres0.K-${K}
BNAME=$outdir/$bname
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

echo "cd $workdir; Rscript $rs $inf $outdir $K $topics $tunemodels $meanmax $cellmin $cellmax" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 6 -m beas -M j.yeung@hubrecht.eu
