#!/bin/sh
# Jake Yeung
# 6-make_LDA_output_matrix_higherK.sh
# Make LDA output matrix as expression matrix.
# Things like log scale should be considered
# 2019-03-07

jmarks="H3K4me1"

thres=0.99
jbin="TRUE"  # I used TRUE to make initial analysis looked good, should try on FALSE to see if it changes?
jcenter="TRUE"

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/lda_to_norm_mat.R"
# fixscript="$wd/scripts/processing/motevo/lib/fix_techrep_sampnames_in_exprs_matrix.R"


[[ ! -d $wd ]] && echo "$wd not found, exiting" && exit 1
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
# [[ ! -e $fixscript ]] && echo "$fixscript not found, exiting" && exit 1

# copy from 5a
prefix="cluster"
suffixCounts="build95_B6_from_traj_CorrPeakFilt.cells_from_bin_analysis"  # MARA analysis now includes exprs matrix with more clels
suffix="build95_B6_from_traj.cells_from_bin_analysis"
jcellmin=0
jcellmax=9999999
highestK="25"
Kvec="$highestK"
jmeanfilt=10

for jmark in $jmarks; do
    inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysis_CorrPeakFilt_1000_${suffix}/lda_outputs.meanfilt_${jmeanfilt}.cellmin_${jcellmin}.cellmax_${jcellmax}.binarize.${jbin}/lda_out_meanfilt.B6-${jmark}.CountThres0.K-${Kvec}.Robj"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_${prefix}_${suffixCounts}/${jmark}/mara_input/count_mats_peaks_norm"
    [[ ! -d $outdir ]] && mkdir -p $outdir  # no need -p because you should have run -p in the previous script making sitecounts?
    outname="hiddenDomains_cellmin_${jcellmin}-cellmax_${jcellmax}-binarize_${jbin}-B6_${jmark}.filt_${thres}.center_${jcenter}_K${highestK}.txt"
    outf=$outdir/$outname

    [[ -e $outf ]] && echo "$outf found, continuing" && continue
    [[ -e $outffixed ]] && echo "$outffixed found, continuing" && continue

    cd $wd
    echo "Getting count matrix"
    if [ $jcenter == "TRUE" ]
    then
        Rscript $rs $inf $outf --thres $thres --center --choose_highest_K
    else
        Rscript $rs $inf $outf --thres $thres --choose_highest_K
    fi
done
