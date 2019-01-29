#!/bin/sh
# Jake Yeung
# 3-get_regions_from_bed.sh
# Get regions from distance to promoter, defined by mm10 liftover swissregulon promoters 
# Merge beds at the end
# 2016-07-22


pyscript="/home/jyeung/projects/tissue_specificity_hogenesch_shellscripts/motevo_dhs_scripts_clean/1-get_gene_region/create_gene_region.py"

# inbed comes from /Home/jyeung/projects/tissue_specificity_hogenesch/make_gene_bed.py
inbed="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/gene_regions/promoter_regions.mm10.bed"
outdir="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/gene_regions"
outbed=$outdir/regions_to_extract.bed
mergedoutbed=$outdir/regions_to_extract.merged.bed

[[ ! -e $pyscript ]] && echo "$pyscript not found, exiting" && exit 1
[[ ! -e $inbed ]] && echo "$inbed not found, exiting" && exit 1

python $pyscript $inbed $outbed --distance 100000

# MERGE
awk '$1!="chrM" {print}' $outbed | sort -k1,1V -k2,2n | bedtools merge -i - > $mergedoutbed

