#!/bin/sh
# Jake Yeung
# 9-deduplicate_bams.sh
#  
# 2019-09-05

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_data/oud3697/raw_demultiplexed"
outmain="$inmain"

for dmain in `ls -d ${inmain}/PZ-*`; do
    echo $dmain
    d=${dmain}/tagged/bwaMapped.sorted.bam
    dout=${dmain}/tagged/bwaMapped.dedup.sorted.bam
    samtools view -b -F 1024 $d > $dout
    samtools index $dout
done
