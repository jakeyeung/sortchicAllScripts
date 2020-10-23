#!/bin/sh
# Jake Yeung
# extend_TSS_bed.sh
#  
# 2020-08-14

inf="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/GRCh38/GCF_000001405.39_GRCh38.p13_genomic.parsed.TSS.txt.merged.chromorenamed.4columns.bed"

outf="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/GRCh38/GCF_000001405.39_GRCh38.p13_genomic.parsed.TSS.txt.merged.chromorenamed.4columns.extended_50000.bed"

dist=25000

awk -v FS="\t" -v OFS="\t" -v dist=$dist '{tss=$2; $2=tss-dist; $3=tss+dist} 1' $inf > $outf
