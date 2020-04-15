# Jake Yeung
# Date of Creation: 2020-04-14
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/check_var_filt_on_umap.R
# Be more stringent with variance filtering, check on UMAP 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(topicmodels)

library(hash)
library(igraph)
library(umap)


# Load data ---------------------------------------------------------------

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


jchromos <- paste("chr", seq(25), sep = "")

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

jprefix <- "/home/jyeung/hub_oudenaarden/jyeung"
jwin <- 50000L
jsuffix <- ""
indir <- file.path(jprefix, paste0("data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZF_AllMerged.winsize_", jwin, jsuffix))

infs <- lapply(jmarks, function(jmark){
  if (jsuffix == ""){
    fname <- paste0("lda_outputs.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.Robj")
  } else if (jsuffix == ".imputevarfilt"){
    fname <- paste0("lda_outputs.counts_table_var_filt.", jmark, ".K-30.binarize.FALSE/ldaOut.counts_table_var_filt.", jmark, ".K-30.Robj")
  }
  inf <- file.path(indir, fname)
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

count.mats <- lapply(jmarks, function(jmark){
  inf.raw <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.winsize_50000.chromo_filtered_by_counts_TAfrac_var.stringent.2020-04-14/count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.rds")
  count.mat <- readRDS(inf.raw)
  # randomly remove 1% of rows
  rows.keep <- sample(x = rownames(count.mat), size = 0.99 * nrow(count.mat), replace = FALSE)
  rows.keep <- gtools::mixedorder(rows.keep)
  count.mat <- count.mat[rows.keep, ]
  return(count.mat)
})

print(infs)

out.lst <- lapply(jmarks, function(jmark){
  inf <- infs[[jmark]]
  load(inf, v=T)
  # fiter out chromosomes
  # chromos.all <- sapply(rownames(count.mat), function(x) strsplit(x, ":")[[1]][[1]])
  # chromos.filt.i <- which(chromos.all %in% jchromos)
  # count.mat.filt <- count.mat[chromos.filt.i, ]
  # load cuts from new obj? 
  count.mat.filt <- count.mats[[jmark]]
  return(list(out.lda = out.lda, count.mat = count.mat.filt))
})


dat.umaps.lst <- lapply(out.lst, function(x){
  tm.result <- posterior(x$out.lda)
  topics.mat <- tm.result$topics
  terms.mat <- tm.result$terms
  
  dat.impute.log <- log2(t(topics.mat %*% terms.mat))
  dat.var <- CalculateVarAll(dat.impute.log, jchromos)
  dat.var.raw <- CalculateVarRaw(x$count.mat, sort.bins = TRUE)
  dat.umap <- DoUmapAndLouvain(topics.mat, jsettings) %>%
    left_join(., dat.var) %>%
    left_join(., dat.var.raw)
})


# Plot  -------------------------------------------------------------------

m.lst <- lapply(dat.umaps.lst, function(dat.umap){
  m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette)
})

print(m.lst)




# Add filtered cells  -----------------------------------------------------

cells.keep <- lapply(jmarks, function(jmark){
  inf.cells <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.winsize_50000.chromo_filtered_by_counts_TAfrac_var.stringent.2020-04-14/count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.colnames")
  dat.cells <- data.frame(mark = jmark, cell = fread(inf.cells, header = FALSE)$V1, stringsAsFactors = FALSE)
  return(dat.cells$cell)
})


# Integrate to dat umap  --------------------------------------------------

dat.umaps.merge.lst <- lapply(jmarks, function(jmark){
  dat.umap <- dat.umaps.lst[[jmark]]
  dat.umap <- dat.umap %>%
    rowwise() %>%
    mutate(is.good = cell %in% cells.keep[[jmark]])
})

m.filt.lst <- lapply(dat.umaps.merge.lst, function(dat.umap){
  m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = is.good)) + 
    geom_point(alpha = 0.25) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~is.good)
})
print(m.filt.lst)

m.cutoffs.lst <- lapply(dat.umaps.merge.lst, function(dat.umap){
  m <- ggplot(dat.umap, aes(x = ncuts.var, y = cell.var.within.sum.norm, color = is.good)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_log10() + scale_y_log10() + facet_wrap(~is.good) + 
    scale_color_manual(values = cbPalette, na.value = "grey85")
})

print(m.cutoffs.lst)

m.cuts.lst <- lapply(dat.umaps.merge.lst, function(dat.umap){
  m <- ggplot(dat.umap, aes(y = ncuts.var, x = ncuts, color = is.good)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_log10() + scale_y_log10() 
})

print(m.cuts.lst)

