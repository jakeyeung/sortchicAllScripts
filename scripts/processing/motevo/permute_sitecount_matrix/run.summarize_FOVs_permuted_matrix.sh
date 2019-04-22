#!/bin/sh
# Jake Yeung
# run.summarize_FOVs_permuted_matrix.sh
#  
# 2019-04-11

inf="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/permute_sitecount_matrix/summarize_FOVs_permuted_matrix.sh"
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
for mark in $marks; do
    # echo $mark
    bash $inf $mark
done
