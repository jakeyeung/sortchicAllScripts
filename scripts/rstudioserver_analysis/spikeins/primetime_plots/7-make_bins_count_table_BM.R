# Jake Yeung
# Date of Creation: 2021-01-08
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/7-make_bins_count_table_BM.R
# Finalize the count tables 50kb for processed at GEO 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)



# Load count tables -------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"
jtype <- "bins"

# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.again"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.final.BMAllmerged.50kb_bins"
dir.create(outdir)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

count.mat.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  if (jmark != "H3K27me3"){
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj"))
  } else {
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.Robj"))
  }
  load(inf.lda, v=T)
  return(count.mat)
})



# Write to output .csv ----------------------------------------------------

for (jmark in jmarks){
  print(jmark)
  fname <- paste0("BoneMarrow_qChIC_cuts_50kb_bins.", jmark, ".csv")
  outf <- file.path(outdir, fname)
  fwrite(as.matrix(count.mat.lst[[jmark]]), file = outf)
}

