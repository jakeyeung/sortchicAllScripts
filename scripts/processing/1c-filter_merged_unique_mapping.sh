#!/bin/sh
# Jake Yeung
# 1c-filter_merged_unique_mapping.sh
# Filter merged bam for highest unique mapping for peak calling purposes 
# 2018-12-18
# Discussion on uniquely mapped reads with BWA
# https://bioinformatics.stackexchange.com/questions/508/obtaining-uniquely-mapped-reads-from-bwa-mem-alignment
# https://bioinformatics.stackexchange.com/questions/2404/problems-to-extract-uniquely-mapping-reads-from-bwa-mem-alignment
# problematic 

inbam="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged/BM.H3K4me1.merged.bam"
outbam="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_rep_merged/BM.H3K4me1.merged.unique.bam"

samtools view -h $inbam | grep -v -e 'XA:Z:' -e 'SA:Z:' | samtools view -b > $outbam


