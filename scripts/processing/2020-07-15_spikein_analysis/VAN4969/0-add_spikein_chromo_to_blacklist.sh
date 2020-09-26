#!/bin/sh
# Jake Yeung
# 0-add_spikein_chromo_to_blacklist.sh
#  
# 2020-08-20

inbed="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/human/ENCFF356LFX.nochr.bed"
[[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1
outbed="/hpc/hub_oudenaarden/jyeung/data/databases/blacklists/human/ENCFF356LFX.nochr.SpikeIns.bed"
[[ -e $outbed ]] && echo "$outbed found, exiting" && exit 1

chromobad="J02459.1"
jstart=0
jend=999999

newline="$chromobad\t$jstart\t$jend"
echo -e $newline > $outbed
cat $inbed >> $outbed
