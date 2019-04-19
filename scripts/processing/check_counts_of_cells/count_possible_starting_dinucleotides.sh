#!/bin/sh
# Jake Yeung
# count_possible_starting_dinucleotides.sh
# Count dinucleotides 
# 2019-04-16


inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_split_by_bc/count_thres-0_build95.withchr"

outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_split_by_bc/count_thres-0_build95.withchr.summaries"

for indir in `ls -d ${inmain}/PZ-BM*`; do
    for b in `ls -d ${indir}/PZ-BM*.bam`; do
        bname=$(basename $b)
        # echo $bname
        outf=$outmain/${bname}.summary.out
        # echo $outf
        # get seequence, get first two letters, sort and count occurrence
        samtools view $b | cut -f10 | cut -b 1-2 | sort | uniq -c > $outf
    done
done

