#!/bin/sh
# Jake Yeung
# 4c-run.summarize_barcodes.sh
# Run summarize barcodes, threshold is a parameter. Set to zero you can remove them later 
# 2018-12-20

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/scripts_analysis/quality_controls/summarize_barcodes.R"

cd $wd

# directories are hardcoded
countthres=0
[[ $countthres != [0-9]* ]] && echo "Must be integer: $countthres" && exit 1

echo "Rscript $rs $countthres"
Rscript $rs $countthres
