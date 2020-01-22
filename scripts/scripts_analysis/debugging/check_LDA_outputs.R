# Jake Yeung
# Date of Creation: 2019-11-26
# File: ~/projects/scchic/scripts/scripts_analysis/debugging/check_LDA_outputs.R
# Check LDA outputs

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(hash)
library(igraph)
library(umap)
library(Matrix)
library(topicmodels)
library(scchicFuncs)

# Load LDA output ---------------------------------------------------------

outdir <- "/Users/yeung/data/scchic/pdfs/stemcell_analysis.redo"

jmark <- "H3K4me3"
jmarks <- c("H3K4me1", "H3K4me3")
jstr <- "Linneg"
jbin <- FALSE

jbins <- c(TRUE, FALSE)

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
for (jmark in jmarks){
  for (jbin in jbins){
    outpdf <- file.path(outdir, paste0("check_LDA_BM.", jmark, ".", jstr, ".", jbin, ".pdf"))
    inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All/lda_outputs.PZ-Bl6-BM-All_Unenriched.", jmark, ".2019-11-22.K-50.binarize.", jbin, "/ldaOut.PZ-Bl6-BM-All_Unenriched.", jmark, ".2019-11-22.K-50.Robj")
    if (!file.exists(inf)){
      print(paste("Skipping", inf))
      next
    }
    
    # inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks/lda_outputs.PZ-Bl6-BM-All_Unenriched.H3K4me1.2019-11-23.K-50.binarize.TRUE/ldaOut.PZ-Bl6-BM-All_Unenriched.H3K4me1.2019-11-23.K-50.Robj"
    # inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks/lda_outputs.PZ-Bl6-BM-All_Unenriched.H3K4me1.2019-11-23.K-50.binarize.FALSE/ldaOut.PZ-Bl6-BM-All_Unenriched.H3K4me1.2019-11-23.K-50.Robj"
    assertthat::assert_that(file.exists(inf))
    load(inf, v=T)
    
    if (class(count.mat.orig) != "logical"){
      print("Using count mat orig")
      count.mat <- count.mat.orig
    }
    print(dim(count.mat))
    
    lsi.out <- RunLSI(as.matrix(count.mat))
    jsettings <- umap.defaults
    jsettings$n_neighbors <- 30
    jsettings$min_dist <- 0.1
    jsettings$random_state <- 123
    umap.out.lsi <- umap(lsi.out$u, config = jsettings)
    dat.umap.long.lsi <- data.frame(cell = rownames(umap.out.lsi$layout), lsi1 = umap.out.lsi$layout[, 1], lsi2 = umap.out.lsi$layout[, 2])
    
    m1 <- ggplot(dat.umap.long.lsi, aes(x = lsi1, y = lsi2)) + geom_point(alpha = 0.5, size = 3) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
      scale_color_manual(values = cbPalette)
    
    if (length(out.lda) > 1){
      out.lda <- out.lda[[1]]
    }
    tm.result <- posterior(out.lda)
    jsettings <- umap.defaults
    jsettings$n_neighbors <- 30
    jsettings$min_dist <- 0.1
    jsettings$random_state <- 123
    umap.out <- umap(tm.result$topics, config = jsettings)
    
    dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
      rowwise() %>%
      mutate(experi = ClipLast(cell, jsep = "-"))
    
    dat.umap.long <- dat.umap.long %>%
      mutate(experi = gsub("-G2", "", experi))
    cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#FF0000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
    
    m2 <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.5, size = 3) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
      scale_color_manual(values = cbPalette)
    
    pdf(outpdf, useDingbats = FALSE)
      print(m1)
      print(m2)
    dev.off()
  }
}





# 
# 
# 
# # do variance
# jchromos <- c(paste("chr", seq(19), sep = ""), "chrX", "chrY")
# jfac <- 10^6
# jpseudo <- 0
# 
# dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms) * jfac + jpseudo)
# cells.var.chromo <- CalculateVarAll(dat.impute.log, jchromos)
# dat.umap.long.var <- left_join(dat.umap.long, cells.var.chromo, by = "cell")
# 
# ggplot(dat.umap.long.var, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_viridis_c(direction = -1) + 
#   facet_wrap(~experi)
