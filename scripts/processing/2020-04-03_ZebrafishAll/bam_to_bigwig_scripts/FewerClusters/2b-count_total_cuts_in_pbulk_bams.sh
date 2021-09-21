#!/bin/sh
# Jake Yeung
# 1b-downsample_bams.sh
#  
# 2020-06-19

# first count reads in all bms

jprefix="imputevarfilt.lessstringent"
mapq=40
ps="/hpc/hub_oudenaarden/jyeung/code_for_analysis/SingleCellMultiOmics.2020-06-03/singlecellmultiomics/bamProcessing/bamFilter.py"
indir="/hpc/hub_oudenaarden/jyeung/data/zebrafish_scchic/bams_tagged_merged_by_marks.split_by_clusters.${jprefix}.mapq_${mapq}.FewerClusters.2020-06-17/DedupR1onlyNoAltHits"

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

for jmark in $jmarks; do
    outf=${indir}/counts_in_bams.${jmark}.txt
    [[ -e $outf ]] && echo "$outf found, exiting" && exit 1
    # mincounts=`cut -d"," $outf -f2 | sort -nk 2 | head -n 1`
    # echo $jmark
    # echo $mincounts
    # continue
    for inbam in `ls -d $indir/PZ-ChIC-ZF_${jmark}*.bam`; do
        samtools idxstats $inbam | cut -f3 | awk -v inbam=$inbam 'BEGIN {total=0} {total += $1} END {print inbam "," total}' >> $outf
    done
done
