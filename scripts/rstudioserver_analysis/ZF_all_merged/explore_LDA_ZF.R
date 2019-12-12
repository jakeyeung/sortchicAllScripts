# Jake Yeung
# Date of Creation: 2019-12-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/explore_LDA_ZF.R
# Plot LDA for different marks ZF

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(topicmodels)

library(scchicFuncs)


# Constants ---------------------------------------------------------------


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")

jbin <- "TRUE"
jsuff <- "AllMerged"
# outdir <- "/home/jyeung/data/from_rstudioserver/scchic/pdfs/"
outdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/LDA_downstream_ZF"

jwin <- 100000L

zscore.cutoff <- "top_500"

# load Chloe data
inf.annot <- paste0("/home/jyeung/hpc/databases/gene_tss/gene_tss.winsize_", jwin, ".species_drerio.nochr.bed")
assertthat::assert_that(file.exists(inf.annot))
inf.WKM <- "/home/jyeung/hpc/scChiC/public_data/Baron_et_al_2019_Zebrafish_WKM/from_macbook/Baron_et_al_pseudobulk_Zebrafish_WKM.rds"



# load tx data ------------------------------------------------------------

# from make_tx_dataset_zebrafish_WKM.R
dat.bulk <- readRDS(inf.WKM)
dat.bulk.mat <- dcast(subset(dat.bulk, select = c(gene, celltype, exprs)), gene ~ celltype, value.var = "exprs")
rownames(dat.bulk.mat) <- dat.bulk.mat$gene; dat.bulk.mat$gene <- NULL

# set up reference data

if (is.character(zscore.cutoff)){
  rnk.cutoff <- as.numeric(strsplit(zscore.cutoff, "_")[[1]][[2]])
  assertthat::assert_that(!is.na(rnk.cutoff))
  dat.bulk.keep <- dat.bulk %>%
    group_by(celltype) %>%
    mutate(rnk = rank(-zscore)) %>%
    filter(rnk < 500)
} else {
  dat.bulk.keep <- dat.bulk %>%
    group_by(gene) %>%
    filter(max(abs(zscore)) > zscore.cutoff)
}

# genes.chic <- sapply(rownames(count.mat), function(x) strsplit(x, ";")[[1]][[2]])
# ref.genes.keep <- intersect(as.character(dat.bulk.keep$gene), genes.chic)
# print(length(ref.genes.keep))


# Load and plot  ----------------------------------------------------------


for (jmark in jmarks){
  print(jmark)
  outf <- paste0("ZF_LDA_output.", jmark, ".binarize_", jbin, ".", jsuff, ".withVar.", Sys.Date(), ".pdf")
  outpdf <- file.path(outdir, outf)
  inf <- paste0("/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_ZFWKM_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.ZF_", 
                jsuff, "_", jmark, 
                ".TAcutoff_0.5.countscutoff_500.binfilt_cellfilt.2019-12-05.K-30.binarize.", jbin, 
                "/ldaOut.ZF_", jsuff, "_", jmark, ".TAcutoff_0.5.countscutoff_500.binfilt_cellfilt.2019-12-05.K-30.Robj")
  assertthat::assert_that(file.exists(inf))
  load(inf, v=T)
  
  # Plot UMAP  --------------------------------------------------------------
  
  tm.result <- posterior(out.lda)
  topics.mat <- tm.result$topics
  
  jsettings <- umap.defaults
  jsettings$n_neighbors <- 30
  jsettings$min_dist <- 0.1
  jsettings$random_state <- 123
  umap.out <- umap(topics.mat, config = jsettings)
  
  dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(experi = ClipLast(cell))
  
  m.umap <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~experi) + ggtitle(jmark)
  
  dat.impute.log <- t(tm.result$topics %*% tm.result$terms)
  jchromos <- seq(25)
  
  dat.var <- CalculateVarAll(dat.impute.log, jchromos)
  
  dat.merge <- left_join(dat.umap.long, dat.var)
  
  m.var <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
    scale_color_viridis_c(direction = -1) + ggtitle(jmark)

  m.var.split <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
    scale_color_viridis_c(direction = -1) + ggtitle(jmark) + facet_wrap(~experi)
  
  
  pdf(outpdf, useDingbats = FALSE)
    print(m.umap)
    print(m.var)
    print(m.var.split)
  dev.off() 
}
