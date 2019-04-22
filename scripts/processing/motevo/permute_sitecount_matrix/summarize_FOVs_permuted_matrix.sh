#!/bin/sh
# Jake Yeung
# summarize_FOVs_permuted_matrix.sh
# Summarize scores  
# 2019-04-11
# 

# jmark="H3K4me1"
jmark=$1
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permutation_analysis_${jmark}/mara_output"
[[ ! -e $indir ]] && echo "$indir not found, exiting" && exit 1

for indir in `ls -d $indir/seed_row_*`; do
    basedir=$(basename $indir)
    jseed=`echo $basedir | cut -d"_" -f3`
    FOVf="$indir/FOV"
    [[ ! -e $FOVf ]] && echo "$FOVf not found, continuing" && continue
    # [[ ! -e $FOVf ]] && echo "$FOVf not found, exiting" && exit 1
    awk -v seed=${jseed} -v mark=${jmark} -F $'\t' 'BEGIN {OFS = FS} {print $1, seed, mark}' ${FOVf}  # print it STDOUT 
    # exit 0
done
