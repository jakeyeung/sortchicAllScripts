#!/bin/sh
# Jake Yeung
# 2-plot_correlation.sh
#  
# 2019-03-21
# https://deeptools.readthedocs.io/en/develop/content/tools/plotCorrelation.html

. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3

jtypes="scatterplot heatmap"
cors="pearson spearman"

suffix="build95"

inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/Lara-Astiaso_2014_${suffix}"
outdir=$inmain

n=0
maxjobs=4

# get cell types and marks from from file name
for f in $(ls -d $inmain/*.npz); do
    # fnames: H3K4me3_CD8_comparison.npz
    # echo $f
    fbase=$(basename $f)
    fbase=${fbase%%.*}
    jmark=$(echo $fbase | cut -d"_" -f1)
    ctype=$(echo $fbase | cut -d"_" -f2)

    for cor in $cors; do
        for jtype in $jtypes; do
            outf=$outdir/${jmark}_${ctype}_comparison.${jtype}.${cor}.pdf
            [[ -e $outf ]] && echo "$outf found, continuing" && continue
            # echo "plotCorrelation --corData $f --corMethod $cor --removeOutliers --skipZeros --whatToPlot $jtype -o $outf --log1p"
            # plotCorrelation --corData $f --corMethod $cor --removeOutliers --skipZeros --whatToPlot $jtype -o $outf --log1p
            # echo $outdir
            # echo ${jmark}_${ctype}_comparison.${jtype}.${cor}.pdf
            echo "plotCorrelation --corData $f --corMethod $cor --removeOutliers --skipZeros --whatToPlot $jtype -o $outf --log1p"
            plotCorrelation --corData $f --corMethod $cor --removeOutliers --skipZeros --whatToPlot $jtype -o $outf --log1p
            if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
                # define maxjobs and n using maxjobsn skeleton
                wait # wait until all have finished (not optimal, but most times good enough)
                echo $n wait
            fi
        done
    done

    
done

