#!/bin/sh
# Jake Yeung
# run_size_estimation.sh
#  
# 2019-05-04

marks="H3K4me3 H3K27me3 H3K9me3"
# mark="H3K4me3"
js="/Users/yeung/software_install_packages/picard.jar"
for mark in $marks; do
		echo $mark
		inf="/Users/yeung/data/scchic/bams_merged/bam_rep_merged/BM_${mark}_merged.bam"
		outf1="/Users/yeung/data/scchic/picarrd_stats/size_metrics.${mark}.txt"
		outf2="/Users/yeung/data/scchic/picarrd_stats/size_metrics.${mark}.pdf"
		[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
		java -jar $js CollectInsertSizeMetrics I=$inf O=$outf1 H=$outf2 M=0.5&
done
