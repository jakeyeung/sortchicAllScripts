# Jake Yeung
# Date of Creation: 2020-01-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/debug/show_var_across_systems_get_objs.R
# description

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)
library(scchicFuncs)

library(parallel)
library(hash)
library(igraph)
library(umap)

# Functions ---------------------------------------------------------------




# Load all intestinal data ------------------------------------------------

jwin <- "50000_25000"
inmain <- "/home/jyeung/hpc/intestinal_scchic/LDA_outputs/topicmodels/ldaAnalysisBins_intestines.2020-01-11.10000_5000"
jbase <- paste0("lda_outputs.mat.Scraped.Unenriched.*.TAcutoff_0.5.countscutoff_500_1000.highcutoff_100000.winsize_", jwin, ".2020-01-11.K-30.binarize.FALSE")
indirs <- list.files(file.path(inmain), pattern = jbase, full.names = TRUE)

infs <- lapply(indirs, function(indir){
  inf <- list.files(indir, pattern = "*.Robj", full.names = TRUE)
  assertthat::assert_that(length(inf) == 1)
  return(inf)
})

jmarks <- sapply(infs, function(x) strsplit(basename(x), split = "\\.")[[1]][[5]])
names(jmarks) <- jmarks
names(infs) <- jmarks

# calculate reads per chromosome per cell, calculate variance, get umap

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

# jmarks <- "k4me3"

out.lst <- mclapply(jmarks, function(jmark){
# out.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- infs[[jmark]]
  load(inf, v=T)  # count.mat, out.lda
  tm.result <- posterior(out.lda)
  dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
  dat.var <- CalculateVarAll(dat.impute.log, jchromos)
  reads.by.chromo <- SumAcrossChromos(count.mat, jchromos, colfunction = mean)
  reads.by.chromo.sum <- SumAcrossChromos(count.mat, jchromos, colfunction = sum)
  reads.by.chromo$label <- jmark
  reads.by.chromo.sum$label <- jmark
  dat.var$label <- jmark
  # put it into one data frame
  ncuts.dat <- data.frame(cell = colnames(count.mat), ncuts = colSums(count.mat), stringsAsFactors = FALSE)
  dat.var <- left_join(dat.var, ncuts.dat)
  return(list(reads.by.chromo = reads.by.chromo, reads.by.chromo.sum = reads.by.chromo.sum, dat.meta = dat.var))
}, mc.cores = length(jmarks))
# })

saveRDS(out.lst, file = "/home/jyeung/data/from_rstudioserver/intestinalchic/var_across_intestinal_marks.rds")

