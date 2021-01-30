# Jake Yeung
# Date of Creation: 2021-01-30
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_H3K9me3/10-check_H3K27me3_features.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)



# Metas  ------------------------------------------------------------------

indir.metas <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage"
dat.metas <- lapply(jmarks, function(jmark){
  fname <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  inf <- file.path(indir.metas, fname)
  fread(inf)
})


# Check pseudbulks of dynamic bins?  -------------------------------------------------------



# load LDA 
outs.all.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  if (jmark != "H3K27me3"){
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj"))
  } else {
    inf.lda <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_H3K27me3.cleaned.varfilt_2.K-30.Robj"))
  }
  load(inf.lda, v=T)
  return(list(out.lda = out.lda, count.mat = count.mat))
})

# jmark.check <- "H3K9me3"
dat.pbulk.lst <- lapply(jmarks, function(jmark.check){
  print(jmark.check)
  cnames.keep.lst <- split(x = dat.metas[[jmark.check]]$cell, f = dat.metas[[jmark.check]]$cluster)
  pbulks.lst <- SumAcrossClusters(outs.all.lst[[jmark.check]]$count.mat, cnames.keep.lst = cnames.keep.lst)
  
  
  mat.pbulk <- bind_rows(pbulks.lst) %>%
    as.data.frame()
  rownames(mat.pbulk) <- rownames(outs.all.lst[[jmark.check]]$count.mat)
  mat.pbulk <- sweep(mat.pbulk, MARGIN = 2, STATS = colSums(mat.pbulk), FUN = "/")
  
  dat.pbulk <- as.matrix(mat.pbulk) %>%
    melt()
  colnames(dat.pbulk) <- c("bin", "ctype", "cuts")
  dat.pbulk$mark <- jmark.check
  return(dat.pbulk)
})




