#!/bin/sh
# Jake Yeung
# 2-plot_correlation.sh
#  
# 2019-03-21
# https://deeptools.readthedocs.io/en/develop/content/tools/plotCorrelation.html


. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3

ctypes="Neu Pro"
marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
jtypes="scatterplot heatmap"
cors="pearson spearman"

suffix="build95_B6"
subdir="log1p_bigwigs"

for ctype in $ctypes; do
    inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/${ctype}_comparison_${suffix}_${subdir}"
    [[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
    outdir=$inmain
    for jmark in $marks; do
        indat="$inmain/${jmark}_${ctype}_comparison.npz"
        [[ ! -e $indat ]] && echo "$indat not found, exiting" && exit 1
        for jtype in $jtypes; do
            for cor in $cors; do
                # echo $jtype $cor
                outf="$outdir/${jmark}_${ctype}_comparison.${jtype}.${cor}.pdf"
                echo "plotCorrelation --corData $indat --corMethod $cor --removeOutliers --skipZeros --whatToPlot $jtype -o $outf"
                plotCorrelation --corData $indat --corMethod $cor --removeOutliers --skipZeros --whatToPlot $jtype -o $outf
            done
        done
    done
done



