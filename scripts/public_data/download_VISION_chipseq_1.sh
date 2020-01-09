#!/bin/sh
# Jake Yeung
# download_VISION_chipseq.sh
#  
# 2019-12-17

inf="http://usevision.org/data/mm10/inputs/"
outdir="/hpc/archive/hub_oudenaarden/public_datasets/BoneMarrow_ChIPseq_Xiang_et_al"
cd $outdir
wget --no-parent -r $inf
