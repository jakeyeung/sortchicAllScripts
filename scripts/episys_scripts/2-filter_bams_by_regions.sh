#!/bin/sh
# Jake Yeung
# 2-filter_bams_by_regions.sh
#  
# 2019-06-25

clstrs="3 5 6"
regions="/Users/yeung/data/scchic/for_episys/regions/regions_to_filter.txt"
[[ ! -e $regions ]] && echo "$regions not found, exiting" && exit 1
outdir="/Users/yeung/data/scchic/for_episys/sorted_bams_filtered"
[[ ! -d $outdir ]] && mkdir $outdir

for clst in $clstrs; do
  inbam="/Users/yeung/data/scchic/for_episys/sorted_bams_build95_B6_stringent/H3K4me3_cluster_${clst}.bam"
  [[ ! -e $inbam ]] && echo "$inbam not found, exiting" && exit 1
  bname=$(basename $inbam)
  bname=${bname%.*}
  outf=$outdir/$bname.filtered.bam
  bedtools intersect -abam $inbam -b $regions > $outf
done

