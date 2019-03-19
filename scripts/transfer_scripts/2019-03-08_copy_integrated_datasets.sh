#!/bin/sh
# Jake Yeung
# 2019-0308_copy_integrated_datasets.sh
#  
# 2019-03-08

# inf="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/integrated_datasets/bm_integrated_H3K4me1_H3K9me3_seurat.rds"
inf="t2:/hpc/hub_oudenaarden/jyeung/data/scChiC/integrated_datasets/integrated_H3K4me1_H3K4me3_seurat.RData"
outdir="/Users/yeung/data/scchic/robjs"

# rsync -avrL --partial --copy-links $inf $outdir
scp $inf $outdir
