#!/bin/sh
# Jake Yeung
# 2019-01-04_copy_binarize_LDA.sh
# This is run on a local computer
# 2019-01-04

# indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/lda_analysis.h3k27me3.broadpeaks/lda_outputs.meanfilt_1.merge_1000_NoM_binarize.cellmin_1000.cellmax_50000"
indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/lda_analysis.h3k27me3.broadpeaks/lda_outputs.meanfilt_1.merge_1000_NoM_cisTopic.cellmin_1000.cellmax_50000"

outdir="/tmp"

scp -r $indir $outdir
