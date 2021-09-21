#!/bin/sh
# Jake Yeung
# 0-index_bams.sh
#  
# 2020-09-13

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/VAN5046_BM/tagged_bams/merged_bams"
cd $indir
for f in `ls -d $indir/*.bam`; do
    samtools index $f&
done
wait
