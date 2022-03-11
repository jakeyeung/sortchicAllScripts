#!/bin/sh
# Jake Yeung
# 1-download_Cusanovich_bams.sh
#  
# 2021-08-25

link1="https://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/BoneMarrow_62016.bam"
link2="https://krishna.gs.washington.edu/content/members/mouse_ATAC_atlas_website/bams/BoneMarrow_62216.bam"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Cusanovich_2018"
cd $outdir

curl -O $link1
curl -O $link2
