#!/bin/sh
# Jake Yeung
# 5b-run_make_sitecount_matrix_from_bed.sh
# Avoid SQL nonesense and just load data to memory 
# 2019-02-04

# jmark="H3K4me1"
jmarks="H3K4me1 H3K4me3"

jsuffix="_mm9_v1"

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/make_sitecount_matrix_from_bed.R"
# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.long.bed"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

jscale=1
jcenter=0
byrow=1

for jmark in $jmarks; do
    echo $jmark
    # inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/${jmark}/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.long.bed"
    inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene${jsuffix}/${jmark}/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.long.bed"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

    outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/${jmark}/mara_input/sitecount_mats${jsuffix}/hiddenDomains_motevo_merged.closest.long.scale_${jscale}.center_${jcenter}.byrow_${byrow}.redo.txt"
    # create the directories beforehand probably
    if [ byrow == 0 ]
    then
    Rscript $rs $inf $outf --scale $jscale --center $jcenter
    else
    Rscript $rs $inf $outf --scale $jscale --center $jcenter --byrow
    fi
done
