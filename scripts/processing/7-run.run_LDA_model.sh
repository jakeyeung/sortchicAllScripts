#!/bin/sh
# Jake Yeung
# 7-run.run_LDA_model.sh
# Run it 
# 2018-12-19

jmem='16G'
jtime='12:00:00'

workdir="/home/hub_oudenaarden/jyeung/projects/scChiC"

cd $workdir

rs="scripts/processing/run_LDA_model.R"

# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats/PZ-BM-H3K4me1.merged.RData"
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats/PZ-BM-H3K4me1.merged.NoCountThres.Robj"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/lda_outputs.meanfilt.NoCountThres"
[[ ! -d $outdir ]] && mkdir $outdir

bname=$(basename $inf)
bname=${bname%%.*}.CountThres0
BNAME=$outdir/$bname
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

# outf=$outdir/$bname.lda_out.RData
# outftune=$outdir/$bname.lda_tune.RData


[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1

K=15
topics="12,15,25,50,100,150"

echo "cd $workdir; Rscript $rs $inf $outdir $K $topics" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 6
# Rscript $rs $inf $outf $outftune
# | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu
