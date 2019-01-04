#!/bin/sh
# Jake Yeung
# 3-run.cisTopic_downstream.sh
# Run
# 2019-01-03

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/cisTopic_downstream.R"
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/lda_analysis.h3k27me3.broadpeaks/lda_outputs.meanfilt_1.merge_1000_NoM_cisTopic.cellmin_1000.cellmax_50000/H3K27me3_NoCountThres_1000_noM_cisTopicOutput.Robj"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/lda_analysis.h3k27me3.broadpeaks/lda_outputs.meanfilt_1.merge_1000_NoM_cisTopic.cellmin_1000.cellmax_50000/downstream"

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1

# need 25 cores here
nthreads=1
jmem='16G'
jtime='12:00:00'
BNAME=$outdir/cistopic_downstream
DBASE=$(dirname "${BNAME}")
[[ ! -d $DBASE ]] && echo "$DBASE not found, exiting" && exit 1
echo "Rscript $rs $inf $outdir" | qsub -l h_rt=${jtime} -l h_vmem=${jmem} -o ${BNAME}.out -e ${BNAME}.err -m beas -M j.yeung@hubrecht.eu -pe threaded $nthreads
