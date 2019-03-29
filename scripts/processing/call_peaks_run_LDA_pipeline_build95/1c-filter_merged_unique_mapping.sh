#!/bin/sh
# Jake Yeung
# 1c-filter_merged_unique_mapping.sh
# Filter merged bam for highest unique mapping for peak calling purposes 
# 2018-12-18
# Discussion on uniquely mapped reads with BWA
# https://bioinformatics.stackexchange.com/questions/508/obtaining-uniquely-mapped-reads-from-bwa-mem-alignment
# https://bioinformatics.stackexchange.com/questions/2404/problems-to-extract-uniquely-mapping-reads-from-bwa-mem-alignment
# problematic 

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged_build95/add_chromo_prefix"

# also filter only 61M, 62M, 63M, 64M, 65M according to /hpc/hub_oudenaarden/avo/scChiC/extract.pl

n=0
maxjobs=4
for b in `ls -d $inmain/*merged.bam`; do
    echo $b
    bname=$(basename $b)
    bname=${bname%.*}
    outbam=$inmain/$bname.unique.bam
    # samtools view -h $b | grep -v -e 'XA:Z:' -e 'SA:Z:' | awk 'BEGIN {OFS="\t"} $5 == "61M" || $5 == "62M" || $5 == "63M" || $5 == "64M" || $5 == "65" {print $0}' | samtools view -b > $outbam&
    samtools view -h $b | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > $outbam&
    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi
done
wait

# b="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged/BM.H3K4me1.merged.bam"
# outbam="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged/BM.H3K4me1.merged.unique.bam"
# 

