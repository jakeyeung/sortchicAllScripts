#!/bin/sh
# Jake Yeung
# 2019-02-15_copy_mara_output.sh
# Copy MARA output 
# 2019-02-15

bin="FALSE"
indir1="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/H3K9me3/mara_output/hiddenDomains_cellmin_100-cellmax_500000-binarize_${bin}-BM_H3K9me3.filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix-"
indir2="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/H3K27me3/mara_output/hiddenDomains_cellmin_100-cellmax_500000-binarize_${bin}-BM_H3K27me3.filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix-"

outdir="/Users/yeung/data/scchic/from_cluster/mara_analysis"

# indirs="$indir2"
indirs="$indir1 $indir2"

for indir in $indirs; do
		echo $indir
		scp -r $indir $outdir
done
