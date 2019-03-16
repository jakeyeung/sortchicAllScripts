#!/bin/sh
# Jake Yeung
# 2019-02-17_copy_LDA_K562_G1_cells.sh
# Analyze G1 cells K562 
# 2019-02-17

# download stuff

inmain="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/LDA_out_K562/ldaAnalysisBins_MetaCell"
outmain="/Users/yeung/data/scchic/from_cluster/LDA_out_K562/."

scp -r $inmain $outmain
