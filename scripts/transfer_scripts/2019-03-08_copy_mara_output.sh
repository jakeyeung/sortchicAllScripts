#!/bin/sh
# Jake Yeung
# 2019-02-17_copy_LDA_K562_G1_cells.sh
# Analyze G1 cells K562 
# 2019-02-17

# download stuff
# inmain="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/LDA_out_K562/ldaAnalysisBins_MetaCell"

mark="H3K4me3"
inmain="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/$mark/mara_output/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_$mark.filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.1000-"
# inmain="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/$mark/mara_output/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_$mark.filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.1000.promoter-"
outmain="/Users/yeung/data/scchic/from_cluster/mara_analysis"

scp -r $inmain $outmain
