#!/bin/sh
# Jake Yeung
# 5-run_make_sitecount_matrix.sh
#  
# 2019-02-04

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/make_sitecount_matrix.R"
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output/motevo_outputs/sql/motevo_merged.closest.long.sqlite3"
wmfile="/hpc/hub_oudenaarden/jyeung/data/databases/WMs/SwissRegulon/mm10_v2_WMs.filt.list"

# outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output/motevo_outputs/sitecount_mats/H3K4me1_sitecount_matrix.GeneID.txt"
# outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/mara_input/count_mats_binned_norm/H3K4me1_sitecount_matrix.GeneID.txt"

outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/mara_input/sitecount_mats/H3K4me1_sitecount_matrix.norm.GeneID.txt"

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"

[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
[[ ! -e $wmfile ]] && echo "$wmfile not found, exiting" && exit 1
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

cd $wd

Rscript $rs $inf $outf --wmlist $wmfile
