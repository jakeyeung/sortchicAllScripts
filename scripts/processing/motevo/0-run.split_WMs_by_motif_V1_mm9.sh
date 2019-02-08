#!/bin/sh
# Jake Yeung
# 0-run.split_WMs_by_motif.sh
# Split WMs into separate motifs
# 2019-01-29

pyscript="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/split_WMs_by_motif.py"
# inf="/hpc/hub_oudenaarden/jyeung/data/databases/WMs/SwissRegulon/mm10_weight_matrices_v2"
inf="/hpc/hub_oudenaarden/jyeung/data/databases/WMs/SwissRegulon/mm9_weight_matrices_v1"
outdir="/hpc/hub_oudenaarden/jyeung/data/databases/WMs/SwissRegulon/mm9_weight_matrices_v1_split"

[[ ! -d $outdir ]] && mkdir $outdir
[[ ! -e $pyscript ]] && echo "$pyscript not found, exiting" && exit 1
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
[[ ! -d $outdir ]] && echo "$outdir not found, exiting" && exit 1

python $pyscript $inf $outdir
