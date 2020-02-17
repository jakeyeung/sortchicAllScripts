# Jake Yeung
# Date of Creation: 2020-02-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged_GLMPCAcorrected_downstream/2-summarize_GLMPCA_UMAPs.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)
library(scchicFuncs)

library(JFuncs)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(hash)
library(igraph)
library(umap)

library(mixtools)

FitMixtureModelLabelCells <- function(topics.mat, topics.keep, jthres = 0.5, show.plots = TRUE, dat.umap.long = NULL){
  # topics.mat: rows are samples, columns are latent variables / topics
  # dat.umap.long requires columns umap1, umap2 and cell name for plotting mixturre fitting output on umap
  assertthat::assert_that(all(topics.keep %in% colnames(topics.mat)))
  
  if (is.null(names(topics.keep))){
    names(topics.keep) <- topics.keep
  }
  # label cells by topic 
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  mm.celltype.lst <- lapply(topics.keep, function(jtopic){
    print(jtopic)
    tvec.raw <- sort(topics.mat[, jtopic])
    # transform
    tvec <- log(tvec.raw / (1 - tvec.raw))
    # xline <- quantile(tvec, )
    mm <- mixtools::normalmixEM(x = tvec, lambda = c(0.9, 0.1), mu = c(-5, -1), sigma = c(2, 1), k = 2)
    (xline <- min(mm$x[which(mm$posterior[, 1] < jthres)]))
    # xline <- -2
    xcells <- names(tvec)[which(tvec > xline)]
    print(paste(length(xcells), "/", length(tvec), "assigned to", jtopic))
    # get topic value for assigned cells
    tvec.raw.filt <- tvec.raw[xcells]
    
    # xline <- max(mm$x[indx.btwn][which(post.filt[, 2] > 0.5)])
    if (show.plots){
      plot(density(tvec), main = paste(jtopic, jthres), xlab = "Log Odds [log(p / (1 - p))]")
      abline(v = xline, col = 'blue')
      plot.mixEM(mm, whichplots = 2, xlab2 = "Log Odds [log(p / (1 - p))]", main2 = paste(jtopic, jthres))
      abline(v = xline, col = 'blue')
    }
    
    cells.keep <- xcells
    if (show.plots){
      m.check <- PlotXYWithColor(dat.umap.long %>% mutate(is.celltype = cell %in% cells.keep), xvar = "umap1", yvar = "umap2", cname = "is.celltype", jtitle = paste(jtopic, jthres), cont.color = FALSE, col.palette = cbPalette)
      print(m.check)
    }
    
    return(list(topic = jtopic, topic.weight = tvec.raw.filt, celltype = xcells, mm = mm, threshold = xline))
  })
}

TidyMixtureModelOutputs <- function(mm.celltype.lst, dedup = TRUE){
  celltypes <- lapply(names(mm.celltype.lst), function(ctype){
    data.frame(cluster = ctype, cell = mm.celltype.lst[[ctype]]$celltype, 
               topic = mm.celltype.lst[[ctype]]$topic, 
               topic.weight = mm.celltype.lst[[ctype]]$topic.weight, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  
  # handle replicates
  if (dedup){
    celltypes <- celltypes %>%
      group_by(cell) %>%
      filter(topic.weight == max(topic.weight))
  }
  return(celltypes)
}

# Constants ---------------------------------------------------------------

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")


# jmark <- "H3K4me1"
# jmark <- "H3K4me1"
jmark <- "H3K9me3"
jexperi <- "Unenriched"
# jexperi <- "AllMerged"

nbins <- "250"
# jcovar.cname <- "ncuts.var"; inf.glm <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2/PZ_", jmark, ".", jexperi, ".KeepBestPlates2.GLMPCA_var_correction.mergebinsize_1000.binskeep_", nbins, ".devmode.2020-02-11.RData")
# jcovar.cname <- "cell.var.within.sum.norm.CenteredAndScaled"
jcovar.cname <- "ncuts.var.log2.CenteredAndScaled"
jpenalty <- 1
# jcovar.cname <- "cell.var.within.sum.norm.CenteredAndScaled"
# inf.glm <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.KeepBestPlates2/PZ_", jmark, ".", jexperi, ".KeepBestPlates2.GLMPCA_var_correction.mergebinsize_1000.binskeep_", nbins, ".covar_", jcovar.cname, ".devmode.2020-02-11.RData")
inf.glm <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.KeepMorePlates/PZ_", jmark, ".", jexperi, ".KeepMorePlates.GLMPCA_var_correction.mergebinsize_1000.binskeep_", nbins, ".covar_", jcovar.cname, ".penalty_", jpenalty, ".winsorize_TRUE.2020-02-04.RData")
assertthat::assert_that(file.exists(inf.glm))

jinf.tss <- "/home/jyeung/hpc/databases/gene_tss/gene_tss_winsize.50000.bed"
assertthat::assert_that(file.exists(jinf.tss))

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Load LDA output ---------------------------------------------------------


inf.lda <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-04.var_filt.UnenrichedAndAllMerged/lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-04.", jexperi, ".K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-04.", jexperi, ".K-30.Robj")
assertthat::assert_that(file.exists(inf.lda))
load(inf.lda, v=T)

tm.result <- posterior(out.lda)
# label  tpic names
colnames(tm.result$topics) <- paste("topic", colnames(tm.result$topics), sep = "")
rownames(tm.result$terms) <- paste("topic", rownames(tm.result$terms), sep = "")

# jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
# dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
# dat.var <- CalculateVarAll(dat.impute.log, jchromos)

# annots.out <- AnnotateBins2(tm.result$terms, top.thres = 0.995, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db")
# annots.out$terms.annot$gene <- sapply(annots.out$terms.annot$termgene, function(x) ifelse(is.na(x), NA, strsplit(x, ";")[[1]][[2]]))
# 
# termgene.dat <- annots.out$terms.annot %>%
#   group_by(term) %>%
#   filter(rnk == min(rnk))

dat.umap.lda <- DoUmapAndLouvain(topics.mat = posterior(out.lda)$topics, jsettings = jsettings) %>%
  rowwise() %>%
  mutate(plate = ClipLast(as.character(cell), jsep = "_"))

# Load GLMPCA correctoin --------------------------------------------------



load(inf.glm, v=T)

tm.result.glm <- list(topics = as.matrix(glm.out$factors), terms = t(as.matrix(glm.out$loadings)))
dat.impute.glm <- t(tm.result.glm$topics %*% tm.result.glm$terms)

# Pload umap before and after  ---------------------------------------------



dat.umap.glm <- DoUmapAndLouvain(topics.mat = glm.out$factors, jsettings = jsettings) %>%
  rowwise() %>%
  mutate(plate = ClipLast(as.character(cell), jsep = "_"))
  
dat.umap.glm$louvain.glm <- dat.umap.glm$louvain; dat.umap.glm$louvain <- NULL

dat.umap.glm <- left_join(dat.umap.glm, subset(dat.umap.lda, select = c(cell, louvain)))


m.lda <- ggplot(dat.umap.lda, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  theme_bw() +  facet_wrap(~plate, ncol = 1) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = cbPalette) + ggtitle("From LDA")
m.glm <- ggplot(dat.umap.glm, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +  
  scale_color_manual(values = cbPalette) + ggtitle("GLMPCA, colored by LDA louvain") + 
  facet_wrap(~plate, ncol = 1)
m.glm.new.louvain <- ggplot(dat.umap.glm, aes(x = umap1, y = umap2, color = louvain.glm)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = cbPalette) + ggtitle("GLMPCA new louvain") + 
  facet_wrap(~plate, ncol = 1)

multiplot(m.lda, m.glm, m.glm.new.louvain, cols = 3)

# 
# 
# # Plot celltypes  ---------------------------------------------------------
# 
# # jtopics <- paste("topic", c(14, 23, 4, 6, 30, 25, 11, 3), sep = "")
# # names(jtopics) <- c("Eryth", "Basophil", "BcellTopic4", "Dendritic", "BcellTopic30", "Neutrophils", "InnateLymph" "Basophil2")
# 
# 
# jtopics <- paste("topic", c(14, 3, 23, 6, 24, 11, 2, 25), sep = "")
# names(jtopics) <- c("Eryth", "Basophil", "Basophil2", "Dendritic", "Bcells", "InnateLymph", "Neutrophils", "Neutrophils2")
# names.final <- c("Eryth", "Basophil", "Basophil", "Dendritic", "Bcells", "InnateLymph", "Neutrophils", "Neutrophils")
# jnames.hash <- hash::hash(names(jtopics), names.final)  # allows clusters to be merged afterwards by having the same final name
# # names(jtopics) <- rep("Eryth", length(jtopics))
# 
# mm.celltype.lst <- FitMixtureModelLabelCells(tm.result$topics, jtopics, jthres = 0.5, show.plots = TRUE, dat.umap.long = dat.umap.lda)
# 
# # plot output
# dat.celltypes <- TidyMixtureModelOutputs(mm.celltype.lst, dedup = TRUE) %>%
#   rowwise() %>%
#   mutate(cluster.orig = cluster,
#          cluster = jnames.hash[[cluster.orig]]) %>%
#   left_join(dat.umap.lda, .)
# 
# ggplot(dat.celltypes, aes(x = umap1, y = umap2, color = cluster)) + geom_point() + 
#   scale_color_manual(values = cbPalette, na.value = "grey85") + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# 
# # Plot labels in GLM corrected dlouvain  ----------------------------------
# 
# 
# 
# 
