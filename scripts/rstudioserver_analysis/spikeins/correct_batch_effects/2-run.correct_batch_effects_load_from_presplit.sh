#!/bin/sh
# Jake Yeung
# 2-run.correct_batch_effects_load_from_presplit.sh
#  
# 2020-12-28

jmark="H3K9me3"
# jmark="H3K27me3"

rs="/home/jyeung/projects/scchic/scripts/rstudioserver_analysis/spikeins/correct_batch_effects/correct_batch_effects_load_from_presplit.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

maindir="/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_batch_correction_output"
inf="${maindir}/imputed_long.${jmark}.from_LDA.binskeep_1000.rds"
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
outf="${maindir}/mat_adj.${jmark}.from_LDA.binskeep_1000.RData"

[[ -e $outf ]] && echo "$outf found, exiting" && exit 1

Rscript $rs $inf $outf -ncores 12

