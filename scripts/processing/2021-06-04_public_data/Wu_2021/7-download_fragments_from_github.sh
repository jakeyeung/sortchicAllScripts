#!/bin/sh
# Jake Yeung
# 7-download_fragments_from_github.sh
#  
# 2021-06-13

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/from_github"

cd $outdir

l1="https://sfurlan.com/cell_lines/sc_fragfiles/K27me3_h1de.fragments.tsv.gz"
l2="https://sfurlan.com/cell_lines/sc_fragfiles/K27me3_stdcells.fragments.tsv.gz"

curl -O $l1
curl -O $l2

