#!/bin/sh
# Jake Yeung
# 0-add_chr_to_bed.sh
#  
# 2021-06-04

inmain="/hpc/hub_oudenaarden/jyeung/data/public_data/ENCODE/bedfiles"

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

for jmark in $jmarks; do
    inbed="${inmain}/ENCODEpeaks.${jmark}.nochr.cleanchr.qval_3.bed"
    [[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1
    outbed="${inmain}/ENCODEpeaks.${jmark}.chr.cleanchr.qval_3.bed"
    [[ -e $outbed ]] && echo "$outbed found, exiting" && exit 1
    sed 's/^/chr/' $inbed > $outbed
done
