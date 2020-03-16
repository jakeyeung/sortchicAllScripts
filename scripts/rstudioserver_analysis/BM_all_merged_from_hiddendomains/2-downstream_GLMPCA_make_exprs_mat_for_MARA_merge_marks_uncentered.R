# Jake Yeung
# Date of Creation: 2020-03-02
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged_from_hiddendomains/2-downstream_GLMPCA_make_exprs_mat_for_MARA_merge_marks.R
# 

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


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

# Load LDA and GLMPCA -----------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
ntopics <- "30"
jbinskeep <- 250

jbase <- paste0("PZ_fromHiddenDomains_", "marks_merged", ".AllMerged.KeepBestPlates2.GLMPCA_novar_correction.binskeep_", jbinskeep, ".ntopics_", ntopics, ".2020-02-11")
# outname <- paste0("countmat_", gsub(".RData$", ".txt", basename(inf.glmpca)))
outname <- paste0("countmat_", jbase, ".txt")
outnamepdf <- paste0("countmat_", jbase, ".pdf")
outdir <- paste0("/home/jyeung/hpc/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/marks_merged/mara_input/count_mats_peaks_unnorm_LDA_merged")
dir.create(outdir, recursive = TRUE)
outf <- file.path(outdir, outname)
outpdf <- file.path(outdir, outnamepdf)

# assertthat::assert_that(!file.exists(outf))
# assertthat::assert_that(!file.exists(outpdf))


dat.imputes <- lapply(jmarks, function(jmark){
 
  print(jmark)
  inf.glmpca <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs_fromhiddendomains.NoVarCorrection.KeepBestPlates2/PZ_fromHiddenDomains_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_novar_correction.mergebinsize_1000.binskeep_", jbinskeep, ".covar_ncuts.var.log2.CenteredAndScaled.penalty_1.winsorize_TRUE.ntopics_", ntopics, ".2020-02-11.RData")
  inf.lda <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_BMAllMerged.2020-02-15.from_hiddendomains/lda_outputs.merged.", jmark, ".minlength_1000.cutoff_analysis.merged.withchr.annotated.K-", ntopics, ".binarize.FALSE/ldaOut.merged.", jmark, ".minlength_1000.cutoff_analysis.merged.withchr.annotated.K-", ntopics, ".Robj")
  inf.lda.bins <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj")
  
  assertthat::assert_that(file.exists(inf.glmpca))
  assertthat::assert_that(file.exists(inf.lda))
  assertthat::assert_that(file.exists(inf.lda.bins))
  
  assertthat::assert_that(endsWith(inf.glmpca, suffix = ".RData"))
  
  # load(inf.glmpca, v=T)
  load(inf.glmpca, v=T)
  # take subset of rows
  rnames.keep <- rownames(glm.inits$Y.filt)
  load(inf.lda, v=T)
  
  tm.result <- posterior(out.lda)
  dat.impute <- log2(t(tm.result$topics %*% tm.result$terms))
  dat.impute.sub <- dat.impute[rnames.keep, ]
  return(data.frame(coord = rownames(dat.impute.sub), dat.impute.sub, stringsAsFactors = FALSE))
})

dat.impute.merged <- Reduce(function(...) merge(..., all = TRUE), dat.imputes)


# Recenter everything and write -------------------------------------------

dat.impute.merged[is.na(dat.impute.merged)] <- 0

mat.merged <- dat.impute.merged[, -1]

geneids <- as.character(dat.impute.merged[, 1])

geneids <- sapply(geneids, function(x) strsplit(x, split = ";")[[1]][[1]])

mat.merged.uncentered.output <- data.frame(Gene.ID = geneids, mat.merged, stringsAsFactors = FALSE)

# write to output
fwrite(mat.merged.uncentered.output, file = outf, row.names = FALSE, sep = "\t")



