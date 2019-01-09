#!/bin/sh
# Jake Yeung
# 4c-run.summarize_barcodes.sh
# Run summarize barcodes, threshold is a parameter. Set to zero you can remove them later 
# 2018-12-20

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/lib/summarize_barcodes.R"

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

cd $wd

# indir, countthres, outdir
cell="K562"  # to get MetaData correct
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/barcode_summaries_K562"
countthres=0
[[ $countthres != [0-9]* ]] && echo "Must be integer: $countthres" && exit 1
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/barcode_summaries_K562/summarized"
chips="H3K27me3,H3K4me3"

Rscript $rs $indir $countthres $outdir $cell $chips
