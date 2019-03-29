#!/bin/sh
# Jake Yeung
# 2-plot_correlation.sh
#  
# 2019-03-21
# https://deeptools.readthedocs.io/en/develop/content/tools/plotCorrelation.html


. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3

ctypes="Multimark"
marks="H3K4me3 H3K27me3 H3K9me3"
markref="H3K4me1"
jtypes="scatterplot heatmap"
cors="pearson spearman"
suffix="log1p_"

n=0
maxjobs=8

for ctype in $ctypes; do
    inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/${ctype}_${suffix}comparison"
    outdir=$inmain
    for jmark in $marks; do
        indat="$inmain/${jmark}_${ctype}_${markref}_vs_${jmark}_${suffix}comparison.npz"
        [[ ! -e $indat ]] && echo "$indat not found, exiting" && exit 1
        for jtype in $jtypes; do
            for cor in $cors; do
                # echo $jtype $cor
                outf="$outdir/${jmark}_${ctype}_${markref}_vs_${jmark}_${suffix}comparison.${jtype}.${cor}.pdf"
                # echo "plotCorrelation --corData $indat --corMethod $cor --removeOutliers --skipZeros --whatToPlot $jtype -o $outf"
                plotCorrelation --corData $indat --corMethod $cor --removeOutliers --skipZeros --whatToPlot $jtype -o $outf&
                if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
                    # define maxjobs and n using maxjobsn skeleton
                    wait # wait until all have finished (not optimal, but most times good enough)
                    echo $n wait
                fi
            done
        done
    done
done
wait
