# Jake Yeung
# Date of Creation: 2022-01-03
# File: ~/projects/scchic/scripts/macbook_analysis_2021/double_incubated_analysis/write_inputs_to_pipeline.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


# Load raw counts ---------------------------------------------------------

indir <- "/Users/yeung/data/scchic/from_cluster_2021/glmpca_objects"

infs <- lapply(jmarks, function(jmark){
  inf.tmp <- file.path(indir, paste0("glmpca.", jmark, ".RData"))
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

glm.inits.lst <- lapply(infs, function(jinf){
  load(jinf, v=T)
  return(glm.inits)
})

lapply(glm.inits.lst, function(jdat) dim(jdat$Y.filt))



# Load raw counts (fromLDA) -----------------------------------------------

infs <- lapply(jmarks, function(jmark){
  if (jmark != "H3K27me3"){
    inf.tmp <- paste0("/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04/lda_outputs.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.binarize.FALSE/ldaOut.count_mat_old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.K-30.Robj")
  } else {
    inf.tmp <- paste0("/Users/yeung/Dropbox/macbookpro_data/data/scchic/from_cluster_2021/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt/lda_outputs.BM_rep2rep3reseq_", jmark, ".cleaned.varfilt_2.K-30.binarize.FALSE/ldaOut.BM_rep2rep3reseq_", jmark, ".cleaned.varfilt_2.K-30.Robj")
  }
  assertthat::assert_that(file.exists(inf.tmp))
  return(inf.tmp)
})

# load(infs[[1]], v=T)

count.mat.lst <- lapply(infs, function(jinf){
  load(jinf, v=T)
  return(count.mat)
})


# Load double incubated ---------------------------------------------------





# Write raw counts rds --------------------------------------------------------




