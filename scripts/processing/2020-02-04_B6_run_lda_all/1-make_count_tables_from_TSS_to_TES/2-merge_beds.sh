#!/bin/sh
# Jake Yeung
# 2-merge_beds.sh
#  
# 2020-08-10

inf1="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTssToTes.chromorenamed.pos.bed.gz"
inf2="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTssToTes.chromorenamed.neg.bed.gz"
[[ ! -e $inf1 ]] && echo "$inf1 not found, exiting" && exit 1
[[ ! -e $inf2 ]] && echo "$inf2 not found, exiting" && exit 1

outf="/hpc/hub_oudenaarden/jyeung/data/databases/refseq/MmRefseqTssToTes.chromorenamed.merged.bed"
outcheck=${outf}.gz
[[ -e $outcheck ]] && echo "$outcheck found, exiting" && exit 1

zcat $inf1 $inf2 | sort -k1,1 -k2,2n > $outf
# zcat $inf2 >> $outf
# sort -k1,1 -k2,2n $outf > $outf
# gzip $outf
