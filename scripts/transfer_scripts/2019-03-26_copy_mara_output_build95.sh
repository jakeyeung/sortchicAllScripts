#!/bin/sh
# Jake Yeung
# 2019-03-26_copy_mara_output_build95.sh
# COpy mara output from build95 
# 2019-03-26

marks="H3K4me3"

outdir="/Users/yeung/data/scchic/from_cluster/mara_analysis_build95/"

for mark in $marks; do
	indir="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95/$mark/mara_output/"
	rsync -avrL --copy-links $indir $outdir
done
