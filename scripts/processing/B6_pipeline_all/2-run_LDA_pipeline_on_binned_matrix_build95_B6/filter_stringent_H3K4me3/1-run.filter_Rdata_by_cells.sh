#!/bin/sh
# Jake Yeung
# 1-run.filter_Rdata_by_cells.sh
#  
# 2019-06-03

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/B6_pipeline_all/2-run_LDA_pipeline_on_binned_matrix_build95_B6/filter_stringent_H3K4me3/filter_Rdata_by_cells.R"
[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA/B6_H3K4me3_pcutoff_0.95_binfilt_cellfilt.2019-05-11.RData"
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
cell="/hpc/hub_oudenaarden/jyeung/data/scChiC/quality_control/good_cells.pcounts.0.15.pfrac.0.05.txt"
[[ ! -e $cell ]] && echo "$cell not found, exiting" && exit 1
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA/stringent_filter/B6_H3K4me3_pcutoff_0.95_binfilt_cellfilt.2019-06-03.RData"
[[ -e $outf ]] && echo "$outf found, exiting" && exit 1

Rscript $rs $inf $cell $outf
