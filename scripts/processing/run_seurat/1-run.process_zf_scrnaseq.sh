#!/bin/sh
# Jake Yeung
# 1-run.process_zf_scrnaseq.sh
#  
# 2019-08-16

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/run_seurat/process_zf_scrnaseq.R"
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/public_data/Zebrafish_WKM/For_Jake/the_massive_complete_zf_dataset.csv.gz"
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/public_data/Zebrafish_WKM/For_Jake/the_massive_complete_zf_dataset.rds"

[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate R3.6; Rscript $rs $inf $outf 
