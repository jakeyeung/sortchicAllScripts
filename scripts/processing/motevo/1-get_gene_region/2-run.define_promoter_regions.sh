#!/bin/sh
# Jake Yeung
# run.define_promoter_regions.sh
# Run /Home/jyeung/projects/sleep_deprivation/motevo/define_promoter_regions.py
# 2016-07-22

bedscript="/Home/jyeung/projects/sleep_deprivation/motevo/define_promoter_regions.py"

[[ ! -e $bedscript ]] && echo "$bedscript not found, exiting" && exit 1

# inf="/archive/epfl/upnae/jyeung/MARA_promoters/mara_promoters_gene_name_association.bed"  # mm9
inf="/archive/epfl/upnae/jyeung/MARA_promoters/mara_promoters_gene_name_association.mm10.bed"
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
outdir="/archive/epfl/upnae/jyeung/sleep_deprivation/motevo_atacseq/gene_regions"
[[ ! -d $outdir ]] && mkdir $outdir
fname="promoter_regions.mm10.bed"
outf=$outdir/$fname

python $bedscript $inf $outf 

