#!/bin/sh
# Jake Yeung
# clean_bed_files.sh
#  
# 2020-10-26

indir="/hpc/hub_oudenaarden/jyeung/data/public_data/ENCODE/bedfiles"

jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

# -log10(qval)
qvalcutoff=3

for jmark in $jmarks; do
    inf="${indir}/ENCODEpeaks.${jmark}.bed"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    outf="${indir}/ENCODEpeaks.${jmark}.nochr.cleanchr.qval_${qvalcutoff}.bed"
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    # sort -k1,1 -k2,2n $inf | sed 's/^chr//g' | awk '$1!~/_/ { print }' > $outf
    awk -v qvalcutoff=$qvalcutoff '$9 > qvalcutoff' $inf | sort -k1,1 -k2,2n  | sed 's/^chr//g' | awk '$1!~/_/ { print }' > $outf
done
