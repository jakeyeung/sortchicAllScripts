#!/bin/sh
# Jake Yeung
# 6-cat_sitecount_matrix.sh
# Run MARA across all marks. Cat sitecount matrices 
# 2020-03-02

jmarks="H3K4me1 H3K4me3 H3K27me3"
# infmain="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged_Peaks_1000"

ndname="sitecount_mats_WithIntercept"
outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/marks_merged/mara_input/${ndname}"
outf="${outdir}/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.addint_1.marks_merged.txt"


[[ -e $outf ]] && echo "$outf found, exiting safetty" && exit 1

for jmark in $jmarks; do
    echo $jmark
    inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/${jmark}/mara_input/${ndname}/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.addint_1.${jmark}.txt"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    # check number of fields in first row
    # awk '{print NF; exit}' $inf
    # continue
    # write to output
    if [ ! -e $outf ]
    then
        echo "Cat to output without removing first row"
    	cat $inf >> $outf  # do not remove first row
    else
        echo "Cat to output while removing first row"
    	cat $inf | sed 1d >> $outf  # remove first row
    fi
done


