#!/bin/sh
# Jake Yeung
# 2019-03-22_copy_LDA_outputs_build95.sh
# Copy output fro build95 and higher count threshold
# 2019-03-22


indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.FALSE.no_filt"
# indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.TRUE.no_filt"
outdir="/Users/yeung/data/scchic/from_cluster/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000/"

# scp -r $indir $outdir

rsync -avrL --copy-links $indir $outdir
