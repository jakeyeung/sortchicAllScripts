#!/bin/sh
# Jake Yeung
# 2019-01-23_copy_LDA_bins.sh
# Copy LDA 
# 2019-01-23

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_MetaCell"
# outdir="/tmp"
outdir="/Users/yeung/data/scchic/from_cluster"

scp -r t2:$indir $outdir
