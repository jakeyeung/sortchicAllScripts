#!/bin/sh
# Jake Yeung
# 1b-handle_overlapping_bins_in_blacklist.sh
#  
# 2020-02-14

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3

inf="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.bed"
infsorted="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.bed"
outfmerged="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/mm10.blacklist.copy.sorted.merged.bed"
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

sort --output $infsorted -k1,1 -k2,2n $inf
bedtools merge -i $infsorted > ${outfmerged}
