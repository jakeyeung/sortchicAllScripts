#!/bin/sh
# Jake Yeung
# 1-project_new_samples_on_LDA.sh
# Project new samples on LDA  
# 2019-06-18

rs="/home/hub_oudenaarden/jyeung/projects/scChiC/scripts/processing/B6_lineage_neg/2-run_LDA/project_new_samples_on_LDA.R"

inlda="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_BinFiltCellFilt2/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.TRUE.no_filt.stringent_filter/lda_out_meanfilt.B6_H3K4me3_pcutoff_0.CountThres0.K-25_30_35_50_mindist_0.4_mindist_processed_lda.Rdata"
inmat="/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-Linneg/PZ-Bl6-BM-Linneg_H3K4me3_binfilt_cellfilt.2019-06-17.RData"

[[ ! -e $rs ]] && echo "$rs not found, exiting" && exit 1
[[ ! -e $inmat ]] && echo "$inmat not found, exiting" && exit 1
[[ ! -e $inlda ]] && echo "$inlda not found, exiting" && exit 1

outf="/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_PZ-Bl6-BM-Linneg/project_from_trained_LDA/lda_out_meanfilt.B6_H3K4me3_pcutoff_0.CountThres0.K-25_30_35_50_mindist_0.4_mindist_processed_lda_vs_PZ-Bl6-BM-Linneg_H3K4me3_binfilt_cellfilt.2019-06-17.RData"

Rscript $rs $inlda $inmat $outf
