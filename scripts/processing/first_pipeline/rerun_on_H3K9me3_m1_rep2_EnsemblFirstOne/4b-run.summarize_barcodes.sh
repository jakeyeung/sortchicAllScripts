#!/bin/sh
# Jake Yeung
# 4b-run.summarize_barcodes.sh
# Run summarize barcodes
# 2019-03-20

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/lib/summarize_barcodes.R"

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

cd $wd

# indir, countthres, outdir
cell="BM"  # to get MetaData correct
# indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/barcode_summaries_K562"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/barcode_summaries"
countthres=0
[[ $countthres != [0-9]* ]] && echo "Must be integer: $countthres" && exit 1
outdir="/home/hub_oudenaarden/jyeung/projects/scChiC/outputs_R/barcode_summaries_Ensembl95"
chips="H3K9me3"

Rscript $rs $indir $countthres $outdir $cell $chips
