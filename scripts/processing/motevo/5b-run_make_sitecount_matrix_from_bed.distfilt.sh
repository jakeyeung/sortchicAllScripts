#!/bin/sh
# Jake Yeung
# 5b-run_make_sitecount_matrix_from_bed.sh
# Avoid SQL nonesense and just load data to memory 
# 2019-02-04

# jmark="H3K4me1"
jmarks="H3K4me1 H3K4me3"
# jmarks="H3K27me3 H3K9me3"

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/motevo/lib/make_sitecount_matrix_from_bed.R"
# inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.long.bed"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

jscale=0
jcenter=0
byrow=0

thres=10000

# fstr="filt.${thres}.promoter"
fstr="filt.${thres}"

for jmark in $jmarks; do
    echo $jmark
    inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/tfbs_output_singlegene/${jmark}/motevo_outputs/bed/merged_bed_closestbed_long/motevo_merged.closest.dist.$fstr.long.bed"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1

    outdir="/hpc/hub_oudenaarden/jyeung/data/scChiC/mara_analysis/${jmark}/mara_input/sitecount_mats"
    [[ ! -d $outdir ]] && mkdir -p $outdir
    outf="${outdir}/hiddenDomains_motevo_merged.closest.long.scale_${jscale}.center_${jcenter}.byrow_${byrow}.bugfix.$fstr.txt"

    # create the directories beforehand probably
    if [ byrow == 0 ]
    then
    Rscript $rs $inf $outf --scale $jscale --center $jcenter
    else
    Rscript $rs $inf $outf --scale $jscale --center $jcenter --byrow
    fi
done
