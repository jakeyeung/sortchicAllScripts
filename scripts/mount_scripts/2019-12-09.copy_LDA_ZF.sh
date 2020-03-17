#!/bin/sh
# Jake Yeung
# 2019-12-09.copy_LDA.sh
# 
# 2019-12-09

inmain="$HOME/hpc/scChiC/raw_demultiplexed/LDA_outputs_all"

# indir1=$inmain/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo
indir2=$inmain/ldaAnalysisBins_ZFWKM_allmarks_mergedtagged_dedupfixed_redo
# indir3=""

outdir="/home/jyeung/data/from_cluster/scchic/LDA_outputs_all"

cp -r $indir2 $outdir
