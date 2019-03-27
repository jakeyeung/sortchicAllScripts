#!/bin/sh
# Jake Yeung
# 0-rename_bigwigs_IdoAmit.sh
# Rename bigwigs 
# 2019-03-26
# Unzip only H3K4me1 and H3K4me3 data

intar="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Lara-Astiaso_2014_Science/GSE60103_RAW.tar"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Lara-Astiaso_2014_Science"

cd $outdir

tar -xf $intar --wildcards "*H3K4me1*"
tar -xf $intar --wildcards "*H3K4me3*"
