#!/bin/sh
# Jake Yeung
# 6-make_LDA_output_matrix.sh
# Make LDA output matrix as expression matrix.
# Things like log scale should be considered
# 2019-02-04

# jmark="H3K4me1"
jmarks="H3K4me1 H3K4me3"

thres=0.99
jbin="FALSE"
jcenter="TRUE"

for jmark in $jmarks; do


    wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
    rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/lda_to_norm_mat.R"
    # inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.FALSE/lda_out_meanfilt.BM-H3K4me1.CountThres0.K-5_15_25.Robj"
    inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisHiddenDomains_1000/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.${jbin}/lda_out_meanfilt.PZ-BM-${jmark}.CountThres0.K-15_20_25_30_35.Robj"

    # outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/count_mats_all/count_mats_binned_norm/cellmin_100-cellmax_500000-binarize_FALSE-BM_H3K4me1.filt_${thres}.txt"
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/${jmark}/mara_input/count_mats_peaks_norm"
    [[ ! -d $outdir ]] && mkdir $outdir
    outname="hiddenDomains_cellmin_100-cellmax_500000-binarize_${jbin}-BM_${jmark}.filt_${thres}.center_${jcenter}.txt"
    outf=$outdir/$outname

    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    # outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/${jmark}/mara_input/count_mats_peaks_norm/hiddenDomains_cellmin_100-cellmax_500000-binarize_FALSE-BM_${jmark}.filt_${thres}.txt"

    cd $wd

    if [ $jcenter == "TRUE" ]
    then
        Rscript $rs $inf $outf --thres $thres --center
    else
        Rscript $rs $inf $outf --thres $thres
    fi
done
