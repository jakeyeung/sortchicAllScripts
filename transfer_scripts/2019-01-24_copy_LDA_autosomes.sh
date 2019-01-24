#!/bin/sh
# Jake Yeung
# 2019-01-24_copy_LDA_autosomes.sh
#  
# 2019-01-24

# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisHiddenDomains_1000/lda_outputs.AutosomesOnly.meanfilt_1.cellmin_100.cellmax_500000.binarize.TRUE"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisHiddenDomains_1000/lda_outputs.AutosomesOnly.meanfilt_1.cellmin_100.cellmax_500000.binarize.FALSE"
outdir="/tmp/ldaAnalysisHiddenDomains_1000"

scp -r t2:$indir $outdir
