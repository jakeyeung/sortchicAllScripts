#!/bin/sh
# Jake Yeung
# 6a-run.fix_techrep_sampnames_in_exprs_matrix.sh
#  
# 2019-03-27

wd="/home/hub_oudenaarden/jyeung/projects/scChiC"

inscript="$wd/scripts/processing/motevo/lib/fix_techrep_sampnames_in_exprs_matrix.R"
[[ ! -e $inscript ]] && echo "$inscript not found, exiting" && exit 1

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
inmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_build95.cells_from_bin_analysis"

# kchoose=50  # maybe bad for H3K4me3??
kchoose="best"

for mark in $marks; do
    indir="${inmain}/${mark}/mara_input/count_mats_peaks_norm"
    fname="hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_${mark}.filt_0.99.center_TRUE_K${kchoose}.txt"
    fnameout="hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_${mark}.filt_0.99.center_TRUE_K${kchoose}_techrepnamefix.txt"

    inf=$indir/$fname
    outf=$indir/$fnameout

    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

    echo "Rscript $inscript $inf $outf"
    cd $wd; Rscript $inscript $inf $outf
done
