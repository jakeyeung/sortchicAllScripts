#!/bin/sh
# Jake Yeung
# 3b-clean_hiddendomains_invalid_records.sh
#  
# 2021-06-06

ps="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/2021-06-04_public_data/remove_invalid_records.py"
[[ ! -e $ps ]] && echo "$ps not found, exiting" && exit 1

inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_Kaya-Okur/H3K27me3/merged_beds/hiddendomains_output/K562_H3K27me3_20181120_allmerged.1000.cutoff/K562_H3K27me3_20181120_allmerged.1000.cutoff_analysis.bed"
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/K562_Kaya-Okur/H3K27me3/merged_beds/hiddendomains_output/K562_H3K27me3_20181120_allmerged.1000.cutoff/K562_H3K27me3_20181120_allmerged.1000.cutoff_analysis.invalid_records_filt.bed"

python $ps $inf $outf
