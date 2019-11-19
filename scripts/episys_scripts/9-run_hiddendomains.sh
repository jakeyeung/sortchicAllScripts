#!/bin/sh
# Jake Yeung
# 9-run_hiddendomains.sh
#  
# 2019-06-26

export PATH="$PATH:/Users/yeung/data/scchic/for_episys/hiddenDomains"

chromsizes="/Users/yeung/data/scchic/for_episys/chromsizes/chromsizes.mm10.filt.txt"
ps="/Users/yeung/data/scchic/for_episys/hiddenDomains/hiddenDomains"

[[ ! -e $chromsizes ]] && echo "$chromsizes not found, exiting" && exit 1

bamdir1="/Users/yeung/data/scchic/for_episys/sorted_bams_filtered"
bamdir2="/Users/yeung/data/scchic/for_episys/sorted_bams_filtered_H3K4me1"

minlength=1000

minpost=0.9
outdir="/Users/yeung/data/scchic/for_episys/hiddenDomains_output"
mapq=30

for bamdir in $bamdir1 $bamdir2; do
	for inf in $(ls -d $bamdir/*.bam); do
		bname=$(basename $inf)
		bname=${bname%%.*}
		bname=${bname}.minpost_${minpost}.mapq_${mapq}
		hiddenDomains -g $chromsizes -b $minlength -t $inf -o $outdir/$bname -p $minpost -q $mapq
	done
done
