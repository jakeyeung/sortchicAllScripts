#!/bin/sh
# Jake Yeung
# 2019-02-08_copy_mara_output.sh
# Copy MARA output 
# 2019-02-08

indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/H3K4me1/mara_output/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_H3K4me1.filt_0.99"

outdir="/Users/yeung/data/scchic/from_cluster/mara_analysis"

scp -r $indir $outdir
