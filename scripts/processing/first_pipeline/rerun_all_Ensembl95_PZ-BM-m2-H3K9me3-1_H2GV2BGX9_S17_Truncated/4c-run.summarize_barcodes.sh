#!/bin/sh
# Jake Yeung
# 4c-run.summarize_barcodes.sh
# Run summarize barcodes, threshold is a parameter. Set to zero you can remove them later 
# 2018-12-20

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/lib/summarize_barcodes.R"

cd $wd

# directories are hardcoded
countthres=0
cell="BM"  # to get MetaData correct
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/barcode_summaries_BM_build95"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/barcode_summaries_BM_build95/summarized"
[[ ! -d $outdir ]] && mkdir $outdir
chips="H3K4me1,H3K4me3,H3K27me3,H3K9me3"

[[ $countthres != [0-9]* ]] && echo "Must be integer: $countthres" && exit 1

echo "Rscript $rs $countthres"
Rscript $rs $indir $countthres $outdir $cell $chips

