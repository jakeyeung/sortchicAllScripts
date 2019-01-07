#!/bin/sh
# Jake Yeung
# 2019-01-04_copy_binarize_LDA.sh
# This is run on a local computer
# 2019-01-04

# indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/lda_analysis.h3k27me3.broadpeaks/lda_outputs.meanfilt_1.merge_1000_NoM_binarize.cellmin_1000.cellmax_50000"
# indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/lda_analysis.h3k27me3.broadpeaks/lda_outputs.meanfilt_1.merge_1000_NoM_cisTopic.cellmin_1000.cellmax_50000"
# jmain="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBroadpeaks_0.5_1000"
# indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBroadpeaks*

indir1="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBroadpeaks_0.5_1000"
indir2="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBroadpeaks_0.3_1000"

for indir in $indir1 $indir2; do
    echo $indir
    bname=$(basename $indir)
    outdir="/tmp/$bname"
    scp -r $indir $outdir
done
