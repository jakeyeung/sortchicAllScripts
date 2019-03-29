#!/bin/sh
# Jake Yeung
# compare_bigwigs.sh
# Compare scChiC and public data 
# 2019-03-21

inf1="/Users/yeung/data/scchic/public_data/Neutrophils/Gong_GenesAndDev_2017/bigwig_mm10/GSM2445206_H3K4me1_Neu_set1.bw"
inf2="/Users/yeung/data/scchic/from_cluster/bigwigs_2019-03-20/bigwigs_2019-03-20/H3K4me1_cluster_3.bw"
outf="/Users/yeung/data/scchic/public_data/bigwig_compares/Neutrophil_Bulk_vs_scChiC.out"

[[ ! -e $inf1 ]] && echo "$inf1 not found, exiting" && exit 1
[[ ! -e $inf2 ]] && echo "$inf2 not found, exiting" && exit 1

bigwigCompare --bigwig1 $inf1 --bigwig2 $inf2 --outFileName $outf
