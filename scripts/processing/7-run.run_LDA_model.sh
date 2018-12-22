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

dist="25000"
meanmax="100"
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats.merge_${dist}/PZ-BM-H3K4me1.merged.NoCountThres.merge_${dist}.Robj"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/lda_outputs.meanfilt_${meanmax}.merge_${dist}"
[[ ! -d $outdir ]] && mkdir $outdir


# outf=$outdir/$bname.lda_out.RData
# outftune=$outdir/$bname.lda_tune.RData


[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1

K=15
topics="10,15,20,50,75,100"
tunemodels="TRUE"

bname=$(basename $inf)
bname=${bname%%.*}.CountThres0.K-${K}
BNAME=$outdir/$bname
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1

echo "cd $workdir; Rscript $rs $inf $outdir $K $topics $tunemodels" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -pe threaded 6 
# Rscript $rs $inf $outf $outftune
# | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu
