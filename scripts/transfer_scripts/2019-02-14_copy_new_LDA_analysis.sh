#!/bin/sh
# Jake Yeung
# 2019-02-14_copy_new_LDA_analysis.sh
# Copy new LDA without filtering binarize TRUE and FALSE 
# 2019-02-14

# indir1="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_10.cellmin_100.cellmax_500000.binarize.TRUE.no_filt"
# indir2="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_10.cellmin_100.cellmax_500000.binarize.FALSE.no_filt"

indir1="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_10.cellmin_100.cellmax_500000.binarize.TRUE"
indir2="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_10.cellmin_100.cellmax_500000.binarize.FALSE"

outdir="/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell"

scp -r t2:$indir1 $outdir/.
scp -r t2:$indir2 $outdir/.
