#!/bin/sh
# Jake Yeung
# 3-make_windows.sh
# Make windows after getting a merged bed 
# 2015-12-03

inbed="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/gene_regions/regions_to_extract.merged.bed"
outbed="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/gene_regions/regions_to_extract.merged.windows500.bed"
bedtools makewindows -b $inbed -w 500 > $outbed
