#!/bin/sh
# Jake Yeung
# copy_matrices_avo_to_jy.sh
# Copy matrices from AvO folder to JY folder for safety 
# 2019-02-15

inmain="/hpc/hub_oudenaarden/avo/scChiC/metacell"
outmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_count_matrix"

cp $inmain/*-100kb.txt $outmain/.
