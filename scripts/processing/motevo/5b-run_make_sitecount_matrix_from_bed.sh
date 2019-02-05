#!/bin/sh
# Jake Yeung
# 5b-run_make_sitecount_matrix_from_bed.sh
# Avoid SQL nonesense and just load data to memory 
# 2019-02-04

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/make_sitecount_matrix_from_bed.R"
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.long.bed"
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.txt"

Rscript $rs $inf $outf
