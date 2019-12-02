#!/bin/sh
# Jake Yeung
# 1-project_new_samples_on_LDA_nobin.sh
# Project new samples on LDA  
# 2019-07-24

jmark="H3K4me1"
rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/B6_lineage_neg/2-run_LDA/project_new_samples_on_LDA_bin.R"

# inlda="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_BinCellFilt2/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.TRUE.no_filt/lda_out_meanfilt.B6_${jmark}_pcutoff_0.CountThres0.K-25_30_40_50_mindist_0.4_mindist_processed_lda.Rdata"
inlda="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/bamlist_for_merging_build95_B6/dat_umap_long_with_louvain.${jmark}.RData"
inmat="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-Linneg/PZ-Bl6-BM-Linneg_${jmark}_binfilt_cellfilt.2019-06-17.RData"
jbin="FALSE"

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
[[ ! -e $inmat ]] && echo "$inmat not found, exiting" && exit 1
[[ ! -e $inlda ]] && echo "$inlda not found, exiting" && exit 1

outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_PZ-Bl6-BM-Linneg/project_from_trained_LDA/bin_${jbin}_PZ-Bl6-BM-Linneg_${jmark}_binfilt_cellfilt.from_louvain.RData"

[[ -e $outf ]] && echo "$outf found, exiting" && exit 1

Rscript $rs $inlda $inmat $outf
