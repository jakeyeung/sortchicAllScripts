#!/bin/sh
# Jake Yeung
# 3-run.downstream_LDA.sh
#  
# 2019-03-22

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/lib/downstream_LDA_100kb.R"

# jmark="H3K4me3"
# jmark="H3K4me1"
# jmark="H3K9me3"
# jmarks="H3K4me1 H3K4me3 H3K9me3"
jmarks="H3K27me3"

n=0
maxjobs=1

for jmark in $jmarks; do
    bin="FALSE"
    Kstr="20_25_30_35"
    suffix=${bin}_${Kstr}
    # inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.${bin}.no_filt/lda_out_meanfilt.BM-${jmark}.CountThres0.K-${Kstr}.Robj"
    # Kstr="15_20_25_30_35"
    inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_10.cellmin_550.cellmax_500000.binarize.${bin}.no_filt/lda_out_meanfilt.BM-${jmark}.CountThres0.K-${Kstr}.Robj"

    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

    plotout="/hpc/hub_oudenaarden/jyeung/data/scChiC/plot_outputs/${jmark}_Build95_meanfilt_10.cellmin_550.cellmax_500000.binarize.${suffix}.pdf"

    echo "Rscript $rs $jmark $inf $plotout"
    cd $wd; Rscript $rs $jmark $inf $plotout&

    if (( $(($((++n)) % $maxjobs)) == 0 )) ; then
        # define maxjobs and n using maxjobsn skeleton
        wait # wait until all have finished (not optimal, but most times good enough)
        echo $n wait
    fi


done
wait
