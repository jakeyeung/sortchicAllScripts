#!/bin/sh
# Jake Yeung
# 1b-reorganize_bed.sh
# Optional: ideally you would incorporate this into step 1 when you cat beds but if you did not you can run this script
# 2016-07-30

inf="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/motevo_merged_bed_long/motevo_merged.closest.long.bed"
outf="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/motevo_merged_bed_long/motevo_merged.closest.long.reorg.bed"
sed 's/\;/\t/g' $inf | sed 's/mm10_//g' | awk -F"\t" -v OFS="\t" ' { t = $5; $5 = $6; $6 = $7; $7 = $8; $8 = t; print; } ' | tr -d $"\r" > $outf
