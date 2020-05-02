#!/bin/sh
# Jake Yeung
# 0-setup_pegasus_data_to_bed.sh
#  
# 2020-04-15

inf="/hpc/hub_oudenaarden/jyeung/data/databases/PEGASUS/danRer7/danRer7_CNEs_PEGASUS.data"
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
outf="/hpc/hub_oudenaarden/jyeung/data/databases/PEGASUS/danRer7/danRer7_CNEs_PEGASUS.forliftover.bed"

awk -F $'\t' 'BEGIN {OFS = FS} NR > 1 {print $1, $2, $3, $4} ' $inf  > $outf
