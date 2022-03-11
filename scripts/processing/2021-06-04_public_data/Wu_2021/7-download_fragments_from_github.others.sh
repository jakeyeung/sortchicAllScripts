#!/bin/sh
# Jake Yeung
# 7-download_fragments_from_github.sh
#  
# 2021-06-13

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/Wu_et_al_2021/from_github"

cd $outdir

# l1="https://sfurlan.com/cell_lines/sc_fragfiles/K27me3_h1de.fragments.tsv.gz"
# l2="https://sfurlan.com/cell_lines/sc_fragfiles/K27me3_stdcells.fragments.tsv.gz"
l1="https://sfurlan.com/pbmc/fragments.tsv.gz"
l2="https://sfurlan.com/gbm/UW7R_raw_data/fragments.tsv.gz"
l3="https://sfurlan.com/gbm/UW7_raw_data/fragments.tsv.gz"

curl $l1 -o "pbmc_fragments.tsv.gz"
curl $l2 -o "UW7R_raw_data_fragments.tsv.gz"
curl $l3 -o "UW7_raw_data__fragments.tsv.gz"

