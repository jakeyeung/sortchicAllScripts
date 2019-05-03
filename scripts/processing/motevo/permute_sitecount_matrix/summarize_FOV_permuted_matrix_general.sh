#!/bin/sh
# Jake Yeung
# 4-summarize_zscores_permuted_matrix_general.sh
# Summarize scores  
# 2019-04-05
# 

jmark=$1
indir=$2
[[ ! -e $indir ]] && echo "$indir not found, exiting" && exit 1

for indir in `ls -d $indir/seed_row_*`; do
    basedir=$(basename $indir)
    jseed=`echo $basedir | cut -d"_" -f3`
    fovf="$indir/FOV"
    [[ ! -e $fovf ]] && echo "$fovf not found, exiting" && exit 
    awk -v seed=${jseed} -v mark=${jmark} -F $'\t' 'BEGIN {OFS = FS} {print $1, seed, mark}' ${fovf}  # print it STDOUT 
    # exit 0
done
