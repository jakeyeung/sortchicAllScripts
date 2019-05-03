#!/bin/sh
# Jake Yeung
# 6-make_LDA_output_matrix_higherK.sh
# Make LDA output matrix as expression matrix.
# Things like log scale should be considered
# 2019-04-16

# jmark="H3K4me1"
# jmarks="H3K4me1 H3K4me3"
# jmarks="H3K27me3 H3K9me3"
jmarks="H3K4me1 H3K4me3"

thres=0.99
jbin="TRUE"  # I used TRUE to make initial analysis looked good, should try on FALSE to see if it changes?
# jbin="FALSE"
jcenter="TRUE"

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
rs="$wd/scripts/processing/motevo/lib/lda_to_norm_mat.R"

# fixscript="$wd/scripts/processing/motevo/run_on_build95/6a-run.fix_techrep_sampnames_in_exprs_matrix.sh"
# fixscript="$wd/scripts/processing/motevo/lib/fix_techrep_sampnames_in_exprs_matrix.R"
# [[ ! -e $fixscript ]] && echo "$fixscript not found, exiting" && exit 1

[[ ! -d $wd ]] && echo "$wd not found, exiting" && exit 1
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

suffix="cluster_build95.withchr.cells_from_bin_analysis"
jdate="2019-04-15"
jcellmin=0
jcellmax=9999999
highestK="best"
jmeanfilt=10
Kvec="25_50"

for jmark in $jmarks; do
    inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisHiddenDomains_1000_${suffix}_${jdate}/lda_outputs.meanfilt_${jmeanfilt}.cellmin_${jcellmin}.cellmax_${jcellmax}.binarize.${jbin}/lda_out_meanfilt.PZ-BM-${jmark}.CountThres0.K-${Kvec}.Robj"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_${suffix}/${jmark}/mara_input/count_mats_peaks_norm"
    [[ ! -d $outdir ]] && mkdir -p $outdir  # no need -p because you should have run -p in the previous script making sitecounts?
    outname="hiddenDomains_cellmin_${jcellmin}-cellmax_${jcellmax}-binarize_${jbin}-BM_${jmark}.filt_${thres}.center_${jcenter}_K${highestK}.txt"
    outf=$outdir/$outname
    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    cd $wd
    echo "Making count matrix"
    if [ $jcenter == "TRUE" ]
    then
        Rscript $rs $inf $outf --thres $thres --center
    else
        Rscript $rs $inf $outf --thres $thres
    fi
done
