#!/bin/sh
# Jake Yeung
# 4-run.downsample_count_matrix.sh
#  
# 2019-06-03

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/B6_pipeline_all/2-run_LDA_pipeline_on_binned_matrix_build95_B6/downsample_count_matrix.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1

counts=700
[[ $counts != [0-9]* ]] && echo "Must be integer: $counts" && exit 1

marks="H3K4me1 H3K4me3 H3K27me3 H3K9me3"
for mark in $marks; do
    inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA/B6_${mark}_pcutoff_0.95_binfilt_cellfilt.2019-05-11.RData"
    outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA/downsampled/B6_${mark}_pcutoff_0.95_binfilt_cellfilt.2019-06-03.downsamp_${counts}.RData"
    [[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
    [[ -e $outf ]] && echo "$outf found, skipping" && continue
    Rscript $rs $inf $outf $counts 
done
