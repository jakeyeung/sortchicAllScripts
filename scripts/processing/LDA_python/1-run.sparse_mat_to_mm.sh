#!/bin/sh
# Jake Yeung
# 1-run.sparse_mat_to_mm.sh
#  
# 2019-07-11

# jmark="H3K36me3"
# jdate="2019-07-11"
# indir="/hpc/hub_oudenaarden/jyeung/data/gastru_scchic/from_macbook/count_mat_binfilt_cellfilt_for_LDA_gastru.namebugfix"
# inf="${indir}/gastru_${jmark}_binfilt_cellfilt.${jdate}.mapq_60.RData"
inf="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA/stringent_filter/B6_H3K4me3_pcutoff_0.95_binfilt_cellfilt.2019-06-03.stringent_filter.RData"
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/mm_files_for_lda/input_data/B6_H3K4me3_pcutoff_0.95_binfilt_cellfilt.2019-06-03.stringent_filter.mm"
rs="/home/hub_oudenaarden/jyeung/projects/scchic_gastru/scripts/processing/2-run_LDA_python/sparse_mat_to_mm.R"

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
[[ ! -e $inf ]] && echo "$inf not found, exiting" && exit 1
[[ -e $outf ]] && echo "$outf found, exiting" && exit 1

Rscript $rs $inf $outf

