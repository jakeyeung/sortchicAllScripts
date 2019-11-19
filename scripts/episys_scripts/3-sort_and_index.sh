#!/bin/sh
# Jake Yeung
# 3-sort_and_index.sh
#  
# 2019-06-25

indir="/Users/yeung/data/scchic/for_episys/sorted_bams_filtered"

for b in `ls -d $indir/*.bam`; do
		samtools index $b
done
