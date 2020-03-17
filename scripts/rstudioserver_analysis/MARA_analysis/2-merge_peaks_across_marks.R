# Jake Yeung
# Date of Creation: 2020-03-03
# File: ~/projects/scchic/scripts/rstudioserver_analysis/MARA_analysis/2-merge_peaks_across_marks.R
# Create a new count table by combining peaks across 3 marks


rm(list=ls())

jstart <- Sys.time()

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)
library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

library(gtools)


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

# Load LDA and GLMPCA -----------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
ntopics <- "30"
jbinskeep <- 250

outdir <- "/home/jyeung/hpc/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.hiddenDomains_output/hd_merged.marks_merged"
outname <- paste0("marks_merged.ntopic_", ntopics, ".binskeep_", jbinskeep, ".txt")
outf <- file.path(outdir, outname)

assertthat::assert_that(!file.exists(outf))

glmpca.imputes <- lapply(jmarks, function(jmark){
  
  print(jmark)
  inf.glmpca <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs_fromhiddendomains.NoVarCorrection.KeepBestPlates2/PZ_fromHiddenDomains_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_novar_correction.mergebinsize_1000.binskeep_", jbinskeep, ".covar_ncuts.var.log2.CenteredAndScaled.penalty_1.winsorize_TRUE.ntopics_", ntopics, ".2020-02-11.RData")
  inf.lda <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_BMAllMerged.2020-02-15.from_hiddendomains/lda_outputs.merged.", jmark, ".minlength_1000.cutoff_analysis.merged.withchr.annotated.K-", ntopics, ".binarize.FALSE/ldaOut.merged.", jmark, ".minlength_1000.cutoff_analysis.merged.withchr.annotated.K-", ntopics, ".Robj")
  inf.lda.bins <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj")
  
  assertthat::assert_that(file.exists(inf.glmpca))
  assertthat::assert_that(file.exists(inf.lda))
  assertthat::assert_that(file.exists(inf.lda.bins))
  
  assertthat::assert_that(endsWith(inf.glmpca, suffix = ".RData"))
  
  load(inf.glmpca, v=T)
  
  tm.result.glm <- list(topics = as.matrix(glm.out$factors), terms = as.matrix(t(glm.out$loadings)))
  dat.impute.glm <- t(tm.result.glm$topics %*% tm.result.glm$terms)
  return(data.frame(coord = rownames(dat.impute.glm), dat.impute.glm, stringsAsFactors = FALSE))
})

peaks.union <- lapply(glmpca.imputes, function(x) rownames(x)) %>%
  unlist() %>%
  gtools::mixedsort() %>%
  unique()

print(length(peaks.union))

coords <- sapply(peaks.union, function(x) strsplit(x, ";")[[1]][[1]], USE.NAMES = FALSE)

# Write bedfile -----------------------------------------------------------

dat.bed <- data.frame(chromo = sapply(coords, function(x) GetChromo(x, add.chr = FALSE)),
                      start = sapply(coords, function(x) GetStart(x)),
                      end = sapply(coords, function(x) GetEnd(x)), 
                      coord = sapply(peaks.union, function(x) strsplit(x, ";")[[1]][[2]]), 
                      stringsAsFactors = FALSE)

fwrite(dat.bed, file = outf, sep = "\t", col.names = FALSE, row.names = FALSE)

