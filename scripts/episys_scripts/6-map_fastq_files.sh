#!/bin/sh
# Jake Yeung
# 5-map_fastq_files.sh
#  
# 2019-06-25

f1="/Users/yeung/data/scchic/for_episys/fastq_full/demultiplexedR1_10000rows.fastq.gz"
f2="/Users/yeung/data/scchic/for_episys/fastq_full/demultiplexedR2_10000rows.fastq.gz"
outf="/Users/yeung/data/scchic/for_episys/bwa_output/bwaMapped.bam"
# ref="/Users/yeung/data/scchic/for_episys/databases/primary_assembly_NOMASK_ERCC92.fa"
ref="/Users/yeung/data/scchic/for_episys/databases/primary_assembly_NOMASK_ERCC92.fa"

bwa mem -t 1 $ref $f1 $f2 | samtools view -Sb - > $outf
