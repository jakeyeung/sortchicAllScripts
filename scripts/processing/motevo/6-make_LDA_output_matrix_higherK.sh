#!/bin/sh
# Jake Yeung
# 6-make_LDA_output_matrix_higherK.sh
# Make LDA output matrix as expression matrix.
# Things like log scale should be considered
# 2019-03-07

# jmark="H3K4me1"
jmarks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"

thres=0.99
jbin="TRUE"  # I used TRUE to make initial analysis looked good, should try on FALSE to see if it changes?
# jbin="FALSE"
jcenter="TRUE"

Esuffix="K500"

for jmark in $jmarks; do
    wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
    rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/lda_to_norm_mat.R"
    inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisHiddenDomains_1000_2019-02-17/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.${jbin}/lda_out_meanfilt.PZ-BM-${jmark}.CountThres0.K-50_100_150_250_500.Robj"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

    # outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats_binned_norm/cellmin_100-cellmax_500000-binarize_FALSE-BM_H3K4me1.filt_${thres}.txt"
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/${jmark}/mara_input/count_mats_peaks_norm"
    [[ ! -d $outdir ]] && mkdir $outdir  # no need -p because you should have run -p in the previous script making sitecounts?
    outname="hiddenDomains_cellmin_100-cellmax_500000-binarize_${jbin}-BM_${jmark}.filt_${thres}.center_${jcenter}_${Esuffix}.txt"
    outf=$outdir/$outname

    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    # outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/${jmark}/mara_input/count_mats_peaks_norm/hiddenDomains_cellmin_100-cellmax_500000-binarize_FALSE-BM_${jmark}.filt_${thres}.txt"
    cd $wd
    if [ $jcenter == "TRUE" ]
    then
        Rscript $rs $inf $outf --thres $thres --center --choose_highest_K
    else
        Rscript $rs $inf $outf --thres $thres --choose_highest_K
    fi
done
