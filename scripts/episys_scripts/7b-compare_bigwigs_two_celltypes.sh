#!/bin/sh
# Jake Yeung
# 7-compare_bigwigs.sh
#  
# 2019-06-25

f1="/Users/yeung/data/scchic/for_episys/sorted_bigwigs_3_clusters/H3K4me3_cluster_3.bw"
f2="/Users/yeung/data/scchic/for_episys/sorted_bigwigs_3_clusters/H3K4me3_cluster_5.bw"
# f3="/Users/yeung/data/scchic/for_episys/sorted_bigwigs_3_clusters/H3K4me3_cluster_6.bw"
# g1="/Users/yeung/data/scchic/for_episys/sorted_bigwigs_3_clusters_H3K4me1/H3K4me1_cluster_11.bw"
g2="/Users/yeung/data/scchic/for_episys/sorted_bigwigs_3_clusters_H3K4me1/H3K4me1_cluster_2.bw"
g3="/Users/yeung/data/scchic/for_episys/sorted_bigwigs_3_clusters_H3K4me1/H3K4me1_cluster_5.bw"

binsize=100000
outf="/Users/yeung/data/scchic/for_episys/bigwig_compare_out_2_celltypes/compare_out.npz"

# multiBigwigSummary bins -b $f1 $f2 $f3 $g1 $g2 $g3 -o $outf --binSize $binsize
multiBigwigSummary bins -b $f1 $f2 $g2 $g3 -o $outf --binSize $binsize
