#!/bin/sh
# Jake Yeung
# 5b-run_make_sitecount_matrix_from_bed.sh
# Avoid SQL nonesense and just load data to memory 
# 2019-02-04

jmark="H3K4me3"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/make_sitecount_matrix_from_bed.R"
# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.long.bed"
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/${jmark}/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.long.bed"
jscale=0
jcenter=0

outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/${jmark}/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_${jscale}.center_${jcenter}.txt"

# create the directories beforehand probably

Rscript $rs $inf $outf
