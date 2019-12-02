#!/bin/sh
# Jake Yeung
# 2-run.run_sctransform.sh
#  
# 2019-08-17

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/run_seurat/run_sctransform.R"
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/public_data/Zebrafish_WKM/For_Jake/the_massive_complete_zf_dataset.rds"
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/public_data/Zebrafish_WKM/For_Jake/the_massive_complete_zf_dataset_seurat.rds"

[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
[[ -e $outf ]] && echo "$outf found, exiting" && exit 1

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf
