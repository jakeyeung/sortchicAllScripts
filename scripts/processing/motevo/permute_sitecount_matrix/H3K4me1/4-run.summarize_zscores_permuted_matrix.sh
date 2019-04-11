#!/bin/sh
# Jake Yeung
# 4-run.summarize_zscores_permuted_matrix.sh
#  
# 2019-04-05

bs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/permute_sitecount_matrix/summarize_zscores_permuted_matrix.sh"
jmark="H3K4me1"

outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis/permutation_analysis_${jmark}/permute_summary"
[[ ! -d $outdir ]] && mkdir $outdir
outf=$outdir/zscore_permute_summary.txt
# echo $outf
bash $bs $jmark > $outf
gzip $outf
