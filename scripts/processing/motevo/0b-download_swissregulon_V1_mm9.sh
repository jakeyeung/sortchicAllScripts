#!/bin/sh
# Jake Yeung
# 0-download_swissregulon_V1_mm9.sh
# Download motif database 
# 2019-01-29

# inf="http://swissregulon.unibas.ch/data/mm10_f5/mm10_weight_matrices_v2"
inf="http://swissregulon.unibas.ch/data/mm9/weight_matrices"

# outdir="/hpc/hub_oudenaarden/jyeung/data/databases/WMs/SwissRegulon"
outdir="/hpc/hub_oudenaarden/jyeung/data/databases/WMs/SwissRegulon"
outf="$outdir/mm9_weight_matrices_v1"

cd $outdir

curl -o $outf $inf
