#!/bin/sh
# Jake Yeung
# 4-run.summarize_zscores_permuted_matrix.sh
#  
# 2019-04-05

# /hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_cluster_build95_CorrPeakFilt.withchr.cells_from_bin_analysis/permutation_H3K4me1/mara_output
bs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/permute_sitecount_matrix/summarize_zscores_permuted_matrix_general.sh"
jmark="H3K4me1"
indir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_cluster_build95_CorrPeakFilt.withchr.cells_from_bin_analysis/permutation_${jmark}/mara_output"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_cluster_build95_CorrPeakFilt.withchr.cells_from_bin_analysis/permutation_${jmark}/permute_summary"
[[ ! -d $indir ]] && echo "$indir not found, exiting" && exit 1
[[ ! -d $outdir ]] && mkdir $outdir

outf=$outdir/zscore_permute_summary.txt
# echo $outf
bash $bs $jmark $indir > $outf
gzip $outf
