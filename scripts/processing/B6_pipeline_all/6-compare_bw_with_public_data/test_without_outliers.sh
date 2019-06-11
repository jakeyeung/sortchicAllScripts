#!/bin/sh
# Jake Yeung
# test_without_outliers.sh
# Run once without outliers to compare why it is so low for Pearson 
# 2019-03-29


. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3
plotCorrelation --corData /hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/MatBcell_comparison_build95/H3K27me3_MatBcell_rep1_comparison.npz --corMethod pearson --whatToPlot scatterplot -o /hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/MatBcell_comparison_build95/H3K27me3_MatBcell_comparison.scatterplot.pearson.KeepOutliers.KeepZeros.pdf --log1p

