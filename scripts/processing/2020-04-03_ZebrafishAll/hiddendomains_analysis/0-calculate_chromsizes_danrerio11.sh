#!/bin/sh
# Jake Yeung
# 0-calculate_chromsizes_danrerio11.sh
#  
# 2020-04-24

fcs="/hpc/hub_oudenaarden/jyeung/software/ucsc_utils/fetchChromSizes"
[[ ! -e $fcs ]] && echo "$fcs not found, exiting" && exit 1

outf="/hpc/hub_oudenaarden/jyeung/data/databases/chromsizes/danRer11.chrom.sizes"

$fcs danRer11 > $outf
