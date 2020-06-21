# Jake Yeung
# Date of Creation: 2020-06-17
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/find_TSS_topics_winsizes_for_heatmap_MouseBM.R
# description


rm(list=ls())

library(hash)
library(ggrastr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(preprocessCore)

library(mixtools)

library(scchicFuncs)
library(JFuncs)

library(topicmodels)


library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)
jorg <- "org.Mm.eg.db"
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")




jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks


# Load DE genes -----------------------------------------------------------

# load this first because it loads a lot of objects, might disuprt things

inf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/de_genes_stringent_objects/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.RData"
load(inf.de, v=T)


# Load LDA r GLMPCA ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden/jyeung/data"

# load GLMPCA from bins 
# jmark <- "H3K4me1"

jexperi <- "AllMerged"
mergesize <- "1000"
nbins <- "1000"
jcovar.cname <- "ncuts.var.log2.CenteredAndScaled"
jpenalty <- 1
ntopics <- 30


out.objs <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.glmpca <- file.path(hubprefix, paste0("scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.good_runs/PZ_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_var_correction.mergebinsize_1000.binskeep_1000.covar_ncuts.var.log2.CenteredAndScaled.penalty_1.winsorize_TRUE.2020-02-11.RData"))
  inf.lda <- file.path(hubprefix, paste0("scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData"))
  inf.lda.bins <- file.path(hubprefix, paste0("scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj"))
  load(inf.glmpca, v=T)
  load(inf.lda, v=T)
  load(inf.lda.bins, v=T)
  
  out <- list(dat.umap.glm.fillNAs = dat.umap.glm.fillNAs, dat.umap.lda = dat.umap.lda, glm.out = glm.out, out.lda = out.lda, count.mat = count.mat)
  return(out)
})

jbins <- out.objs$H3K4me1$out.lda@terms

# get imputed mats

dat.imputes.lst <- lapply(out.objs, function(x){
  tm.result <- topicmodels::posterior(x$out.lda)
  dat.impute <- log2(t(tm.result$topics %*% tm.result$terms) * 10^6)
  return(dat.impute)
})

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.vars.lst <- lapply(dat.imputes.lst, function(dat.impute.log){
  CalculateVarAll(dat.impute.log, jchromos)
  
})

