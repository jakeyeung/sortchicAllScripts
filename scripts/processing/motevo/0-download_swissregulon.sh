#!/bin/sh
# Jake Yeung
# 0-download_swissregulon.sh
# Download motif database 
# 2019-01-29

inf="http://swissregulon.unibas.ch/data/mm10_f5/mm10_weight_matrices_v2"

outdir="/hpc/hub_oudenaarden/jyeung/data/databases/WMs/SwissRegulon"

cd $outdir

curl -O $inf
