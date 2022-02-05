#!/bin/sh
# Jake Yeung
# 1-download_meta.sh
#  
# 2021-08-25

# link="http://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/metadata/cell_metadata.txt"
link="https://krishna.gs.washington.edu/content/members/ajh24/mouse_atlas_data_release/metadata/cell_metadata.txt"

indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018/metadata"

cd $indir
curl -O $link
