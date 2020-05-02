#!/bin/sh
# Jake Yeung
# 4-remove_chr_from_bed.sh
# bam files have no chr so remove it from bed
# 2020-04-15

dist=10000
inbed="/hpc/hub_oudenaarden/jyeung/data/databases/PEGASUS/danRer7/danRer11_CNEs_PEGASUS.forliftover.${dist}.bed"
outbed="/hpc/hub_oudenaarden/jyeung/data/databases/PEGASUS/danRer7/danRer11_CNEs_PEGASUS.forliftover.${dist}.nochr.bed"

# awk -v dist=$dist -F $'\t' 'BEGIN {OFS = FS} {a=int(($2+$3)/2); print $1, a - dist/2, a + dist/2, $4  } ' $inf > $outf
# awk -F"\t" 'BEGIN {OFS = FS} {print substr($2$1"|"substr($2,2)"|"$3"|"$4}' $inbed > $outbed

# remove chr, then check first column is digit
# awk -F, '$1 ~ /^[[:digit:]]+$/'

sed 's/^chr//g' $inbed | awk -F"\t" '$1 ~ /^[[:digit:]]+$/' > $outbed
