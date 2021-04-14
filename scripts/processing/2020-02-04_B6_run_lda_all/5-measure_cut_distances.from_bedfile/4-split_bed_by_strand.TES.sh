#!/bin/sh
# Jake Yeung
# 4-split_bed_by_strand.sh
#  
# 2020-08-06

inf="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTes.chromorenamed.bed.gz"

outf1="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTes.chromorenamed.pos.bed.gz"
outf2="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTes.chromorenamed.neg.bed.gz"

zcat $inf | awk '$4 == "+" {print $0}' | gzip  > $outf1
zcat $inf | awk '$4 == "-" {print $0}' | gzip  > $outf2
