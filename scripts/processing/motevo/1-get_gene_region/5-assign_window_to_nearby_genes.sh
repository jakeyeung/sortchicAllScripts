#!/bin/sh
# Jake Yeung
# 5-assign_window_to_nearby_genes.sh
# Assign each window to genes based on distance (multiple genes) 
# 2016-04-04

dist=50000
[[ $dist != [0-9]* ]] && echo "Must be integer: $dist" && exit 1

assignscript="/Home/jyeung/projects/tissue_specificity_hogenesch_shellscripts/merged_dhs/expand_gene_collect_peaks.py"
regionsbed="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/gene_regions/regions_to_extract.merged.windows500.bed"
refbed="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/gene_regions/promoter_regions.mm10.bed"  # from swissregulon promoters
outf="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/gene_regions/regions_to_extract.merged.windows500.distfilt.$dist.bed"

[[ ! -e $assignscript ]] && echo "$assignscript not found, exiting" && exit 1
[[ ! -e $regionsbed ]] && echo "$regionsbed not found, exiting" && exit 1
[[ ! -e $refbed ]] && echo "$refbed not found, exiting" && exit 1

python $assignscript $refbed $regionsbed $dist $outf
