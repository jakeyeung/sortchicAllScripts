#!/bin/sh
# Jake Yeung
# 2-plot_correlation.sh
#  
# 2019-03-21
# https://deeptools.readthedocs.io/en/develop/content/tools/plotCorrelation.html


. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3

# ctypes="Neu Pro"
ctypes="MatBcell HSC ProB"
marks="H3K4me1 H3K4me3 H3K27me3"
jtypes="scatterplot heatmap"
cors="pearson spearman"

for ctype in $ctypes; do
    inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/${ctype}_comparison"
    outdir=$inmain
    for jmark in $marks; do
        if [[ $jmark == "H3K27me3" && $ctype == "HSC" ]]; then
            echo "SKipping $jmark $ctype"
            continue
        fi
        indat="$inmain/${jmark}_${ctype}_${ctype}_comparison.npz"
        [[ ! -e $indat ]] && echo "$indat not found, exiting" && exit 1
        for jtype in $jtypes; do
            for cor in $cors; do
                # echo $jtype $cor
                outf="$outdir/${jmark}_${ctype}_comparison.${jtype}.${cor}.pdf"
                [[ -e $outf ]] && echo "$outf  found, continuing" && continue
                echo "plotCorrelation --corData $indat --corMethod $cor --removeOutliers --skipZeros --whatToPlot $jtype -o $outf --log1p"
                plotCorrelation --corData $indat --corMethod $cor --removeOutliers --skipZeros --whatToPlot $jtype -o $outf --log1p
            done
        done
    done
done



