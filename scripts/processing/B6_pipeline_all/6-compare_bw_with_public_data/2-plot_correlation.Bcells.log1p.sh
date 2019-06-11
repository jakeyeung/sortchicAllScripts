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

suffix="build95_B6"
subdir="log1p_bigwigs"

n=0
maxjobs=8

for ctype in $ctypes; do
    inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/public_data/comparisons_with_pseudobulk/${ctype}_comparison_${suffix}_${subdir}"
    [[ ! -d $inmain ]] && echo "$inmain not found, exiting" && exit 1
    outdir=$inmain
    for jmark in $marks; do
        if [[ $jmark == "H3K27me3" && $ctype == "HSC" ]]; then
            echo "SKipping $jmark $ctype"
            continue
        fi
        for indat in `ls -d $inmain/${jmark}_${ctype}_*_comparison.npz`; do
            for jtype in $jtypes; do
                for cor in $cors; do
                    # echo $jtype $cor
                    outf="$outdir/${jmark}_${ctype}_comparison.${jtype}.${cor}.pdf"
                    # [[ -e $outf ]] && echo "$outf  found, continuing" && continue
                    echo "plotCorrelation --corData $indat --corMethod $cor --removeOutliers --skipZeros --whatToPlot $jtype -o $outf"
                    plotCorrelation --corData $indat --corMethod $cor --removeOutliers --skipZeros --whatToPlot $jtype -o $outf&
                    # exit 0
                    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
                    	# define maxjobs and n using maxjobsn skeleton
                        wait # wait until all have finished (not optimal, but most times good enough)
                        echo $n wait
                    fi

                    
                done
            done
        done
        # indat="$inmain/${jmark}_${ctype}_${ctype}_comparison.npz"
        # [[ ! -e $indat ]] && echo "$indat not found, exiting" && exit 1
    done
done
wait


