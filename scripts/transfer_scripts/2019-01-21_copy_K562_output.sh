#!/bin/sh
# Jake Yeung
# 2019-01-21_copy_K562_output.sh
# Copy it  
# 2019-01-21

indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/LDA_out_K562/ldaAnalysishiddenDomains_1000_round2/lda_outputs.hiddenDomains.meanfilt_1.cellmin_100.cellmax_500000.binarize.FALSE"
outdir="/tmp/count_mat_K562_round2_LDA_output"

scp -r $indir $outdir
