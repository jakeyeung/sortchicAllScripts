#!/bin/sh
# Jake Yeung
# 0-download_Cusanovich_rds.sh
#  
# 2021-08-26

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/processed_data"

# link="http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/matrices/atac_matrix.binary.qc_filtered.rds"
link="https://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/matrices/atac_matrix.binary.qc_filtered.rds"

cd $indir
curl -O $link

