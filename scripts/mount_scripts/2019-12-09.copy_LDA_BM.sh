#!/bin/sh
# Jake Yeung
# 2019-12-09.copy_LDA_BM.sh
#  
# 2019-12-09

inmain="$HOME/hpc/scChiC/raw_demultiplexed/LDA_outputs_all"
# indir2=$inmain/ldaAnalysisBins_BM_allmarks_mergedtagged_dedupfixed_redo
indir2=$inmain/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo
[[ ! -d $indir2 ]] && echo "$indir2 not found, exiting" && exit 1

outdir="/home/jyeung/data/from_cluster/scchic/LDA_outputs_all"

cp -r $indir2 $outdir
