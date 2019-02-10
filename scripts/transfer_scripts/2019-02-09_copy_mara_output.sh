#!/bin/sh
# Jake Yeung
# 2019-02-08_copy_mara_output.sh
# Copy MARA output 
# 2019-02-08

# indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/H3K4me1/mara_output/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_H3K4me1.filt_0.99"
# indir1="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/H3K4me3/mara_output/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_H3K4me3.filt_0.99-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0-"
indir1="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/H3K4me3/mara_output/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_H3K4me3.filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0-"
indir2="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/H3K4me1/mara_output/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_H3K4me1.filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0-"

outdir="/Users/yeung/data/scchic/from_cluster/mara_analysis"

# indirs="$indir2"
indirs="$indir1 $indir2"

for indir in $indirs; do
		echo $indir
		scp -r $indir $outdir
done
