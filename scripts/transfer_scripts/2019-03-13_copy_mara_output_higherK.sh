#!/bin/sh
# Jake Yeung
# 2019-03-13_copy_mara_output_higherK.sh
# Copy more mara outputs 
# 2019-03-13

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
suffixE="K500"

# mark="H3K4me1"

for mark in $marks; do
		echo $mark
		indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/${mark}/mara_output/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_${mark}.filt_0.99.center_TRUE_${suffixE}-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0--${suffixE}"
		outmain="/Users/yeung/data/scchic/from_cluster/mara_analysis"
		scp -r $indir $outmain
done
