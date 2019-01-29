#!/bin/sh
# Jake Yeung
# run.downstream_LDA_100kb.sh
# Run Rscript from the correct directory to get the relative paths correct 
# 2019-01-29

inscript="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/scripts_analysis/primetime_from_server/downstream_LDA_100kb.R"
wd="/home/hub_oudenaarden/jyeung/projects/scChiC"

cd $wd

Rscript $inscript
