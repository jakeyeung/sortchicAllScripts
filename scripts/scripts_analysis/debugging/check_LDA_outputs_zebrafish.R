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



# Constants ---------------------------------------------------------------



cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")


# Load LDA output ---------------------------------------------------------

# TSS
inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZF_All_TSS/lda_outputs.PZ-ChIC-ZFWKM-H3K4me3.winsize_50000.merged.K-50.binarize.FALSE/ldaOut.PZ-ChIC-ZFWKM-H3K4me3.winsize_50000.merged.K-50.Robj"

# bins
jmark <- "H3K4me3"
jstr <- "UnenrichedXStemCells"
jbin <- "TRUE"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jstrs <- c("Merged", "UnenrichedXStemCells", "Unenriched", "StemCells")
jbins <- c(TRUE, FALSE)
outdir <- "/Users/yeung/data/scchic/pdfs/zebrafish.redo/try2"

for (jmark in jmarks){
  for (jstr in jstrs){
    for (jbin in jbins){
      inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZF_All/lda_outputs.PZ-ZF-All_", jstr, ".", jmark, ".2019-11-22.K-50.binarize.", jbin, "/ldaOut.PZ-ZF-All_", jstr, ".", jmark, ".2019-11-22.K-50.Robj")
      if (!file.exists(inf)){
        print(paste("Skipping inf", inf))
        next
      }
      outpdf <- file.path(outdir, paste("ZF", jmark, jstr, jbin, "2019-11-26.pdf", sep = "_"))
      
      load(inf, v=T)
      print(dim(count.mat))
      
      if (class(count.mat.orig) == "logical"){
        lsi.out <- RunLSI(as.matrix(count.mat))
      } else {
        lsi.out <- RunLSI(as.matrix(count.mat.orig))
      }
      jsettings <- umap.defaults
      jsettings$n_neighbors <- 30
      jsettings$min_dist <- 0.1
      jsettings$random_state <- 123
      umap.out.lsi <- umap(lsi.out$u, config = jsettings)
      
      
      dat.umap.long.lsi <- data.frame(cell = rownames(umap.out.lsi$layout), pc1 = umap.out.lsi$layout[, 1], pc2 = umap.out.lsi$layout[, 2], stringsAsFactors = FALSE) %>%
        rowwise() %>%
        mutate(is.stem = grepl("CD41plus", cell),
               experi = ClipLast(cell))
      
      m.lsi1 <- ggplot(dat.umap.long.lsi, aes(x = pc1, y = pc2, color = is.stem)) + geom_point(alpha = 0.5, size = 3) + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
        scale_color_manual(values = cbPalette)
      
      m.lsi2 <- ggplot(dat.umap.long.lsi, aes(x = pc1, y = pc2, color = is.stem)) + geom_point(alpha = 0.5, size = 3) + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
        scale_color_manual(values = cbPalette) + 
        facet_wrap(~experi)
      
      if (length(out.lda) > 1){
        out.lda <- out.lda[[1]]
      }
      
      tm.result <- posterior(out.lda)
      jsettings <- umap.defaults
      jsettings$n_neighbors <- 30
      jsettings$min_dist <- 0.1
      jsettings$random_state <- 123
      umap.out <- umap(tm.result$topics, config = jsettings)
      
      m.lda1 <- dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
        rowwise() %>%
        mutate(experi = ClipLast(cell, jsep = "-"))
      
      dat.umap.long <- dat.umap.long %>%
        mutate(experi = gsub("-G2", "", experi))
      cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#FF0000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
      
      m.lda2 <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.5, size = 3) +
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
        scale_color_manual(values = cbPalette)
      
      
      # do variance
      jchromos <- c(paste("chr", seq(19), sep = ""), "chrX", "chrY")
      jfac <- 10^6
      jpseudo <- 0
      
      dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms) * jfac + jpseudo)
      cells.var.chromo <- CalculateVarAll(dat.impute.log, jchromos)
      dat.umap.long.var <- left_join(dat.umap.long, cells.var.chromo, by = "cell")
      
      m.var <- ggplot(dat.umap.long.var, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        scale_color_viridis_c(direction = -1) + 
        facet_wrap(~experi)
      
      pdf(outpdf, useDingbats = FALSE)
      print(m.lsi1)
      print(m.lsi2)
      print(m.lda1)
      print(m.lda2)
      print(m.var)
      dev.off()
    }
  }
}




