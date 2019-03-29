#!/bin/sh
# Jake Yeung
# 2-plot_correlation.sh
#  
# 2019-03-21
# https://deeptools.readthedocs.io/en/develop/content/tools/plotCorrelation.html


. /hpc/hub_oudenaarden/jyeung/software/anaconda3/etc/profile.d/conda.sh; conda activate py3

ctypes="Erythrobl Megakaryo"
marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
reps="rep1 rep2"
jtypes="scatterplot heatmap"
cors="pearson spearman"

n=0
maxjobs=8

for ctype in $ctypes; do
    inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/${ctype}_comparison"
    outdir=$inmain
    for jmark in $marks; do
        if [[ $ctype == "Megakaryo" && $jmark == "H3K9me3" ]]; then
        	echo "Skip $ctype and $jmark"
            continue
        fi
        for rep in $reps; do
            indat="$inmain/${jmark}_${ctype}_${rep}_comparison.npz"
            [[ ! -e $indat ]] && echo "$indat not found, exiting" && exit 1
            for jtype in $jtypes; do
                for cor in $cors; do
                    # echo $jtype $cor
                    outf="$outdir/${jmark}_${ctype}_${rep}_comparison.${jtype}.${cor}.pdf"
                    [[ -e $outf ]] && echo "$outf found, continuing" && continue
                    echo "plotCorrelation --corData $indat --corMethod $cor --removeOutliers --skipZeros --whatToPlot $jtype -o $outf --log1p"
                    plotCorrelation --corData $indat --corMethod $cor --removeOutliers --skipZeros --whatToPlot $jtype -o $outf --log1p
                    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
                    	# define maxjobs and n using maxjobsn skeleton
                        wait # wait until all have finished (not optimal, but most times good enough)
                        echo $n wait
                    fi
                done
            done
        done
    done
done
wait
