#!/bin/sh
# Jake Yeung
# make_binned_matrix.sh
#  
# 2019-08-01

cd "/Users/yeung/data/scchic/bigwigs/antibody_tests"

multiBigwigSummary bins -b PZ-K4me1-85-1-2_100kb.bw PZ-K4me1-85-1-4_100kb.bw PZ-K4me1-96-2-2_100kb.bw PZ-K4me1-96-2-4_100kb.bw -o PZ_K4me1_four_results_10-kb.npz --outRawCounts PZ_K4me1_four_results_100kb.txt --blackListFileName /Users/yeung/data/databases/blacklist_regions/hg19-blacklist.v2.bed --binSize 100000
