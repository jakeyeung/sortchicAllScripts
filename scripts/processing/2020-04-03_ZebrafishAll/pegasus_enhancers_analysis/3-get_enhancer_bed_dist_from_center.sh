#!/bin/sh
# Jake Yeung
# 3-get_enhancer_bed_dist_from_center.sh
#  
# 2020-04-15

# dist=10000
dist=20000
# halfdist=5000
# dist=$(( $halfdist + $halfdist ))
inf="/hpc/hub_oudenaarden/jyeung/data/databases/PEGASUS/danRer7/danRer11_CNEs_PEGASUS.forliftover.bed"
outf="/hpc/hub_oudenaarden/jyeung/data/databases/PEGASUS/danRer7/danRer11_CNEs_PEGASUS.forliftover.${dist}.bed"

awk -v dist=$dist -F $'\t' 'BEGIN {OFS = FS} {a=int(($2+$3)/2); min = a<5000 ? 0 : a - dist/2; max = a + dist/2; print $1, min, max, $4} ' $inf > $outf
