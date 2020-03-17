# Jake Yeung
# Date of Creation: 2020-02-13
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged_GLMPCAcorrected_downstream/3-explore_simplex_of_celltypes.R
# Explore simplelx of celltypes

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(hash)
library(igraph)
library(umap)

library(ggrepel)

GetCondFromSamp <- function(samp, mark = "H3K4me1"){
  if (mark == "H3K4me1"){
    wt.plates <- c("B6-13W1-BM-H3K4me1", "PZ-ChIC-Bl6-BM-H3K4me1-Index")
    linneg.plates <- c("PZ-ChIC-Bl6-BM-lin-H3K4me1-")
    hsc.plates <- c("-stem-cells-")
  } else if (mark == "H3K4me3"){
    wt.plates <- c("B6-13W1-BM-H3K4me3")
    linneg.plates <- c("-Linneg-")
    hsc.plates <- c("-B6BMSC-") 
  } else if (mark == "H3K27me3") {
    wt.plates <- c("B6-13W1-BM-")
    linneg.plates <- c("-Linneg-")
    hsc.plates <- c("-B6BMSC-")  
  } else if (mark == "H3K9me3") {
    wt.plates <- c("B6-13W1-BM-")
    linneg.plates <- c("-Linneg-")
    hsc.plates <- c("BMSC-")  
  } else {
    print(paste(mark, "not yet coded"))
  }
  wt.plates.grep <- paste(wt.plates, collapse = "|")
  linneg.plates.grep <- paste(linneg.plates, collapse = "|")
  hsc.plates.grep <- paste(hsc.plates, collapse = "|")
  is.wt <- grepl(wt.plates.grep, samp)
  is.linneg <- grepl(linneg.plates.grep, samp)
  is.hsc <- grepl(hsc.plates.grep, samp)
  bool.vec <- c(is.wt, is.linneg, is.hsc)
  assertthat::assert_that(sum(bool.vec) == 1)
  
  indx <- which.max(bool.vec)
  if (indx == 1){
    cond <- "Unenriched"
  } else if (indx == 2){
    cond <- "Linneg" 
  } else if (indx == 3){
    cond <- "StemCell"
  } else {
    stop("must be 1, 2, or 3. Found: ", indx)
  }
}


CalculateCelltypeFractions <- function(dat.umap, jfrac = 0.8, jseed = 0){
  set.seed(jseed)
  
  dat.total <- dat.umap %>%
    group_by(cond) %>%
    summarise(ncell = length(cell))
  (jsize <- ceiling(min(dat.total$ncell) * 0.8))
  
  jsum <- dat.umap %>%
    group_by(cond) %>%
    sample_n(tbl = ., size = jsize, replace = FALSE) %>%
    group_by(cond, cluster) %>%
    summarise(ncell = length(cell)) %>%
    group_by(cluster) %>%
    mutate(frac = ncell / sum(ncell))
  
  jsum.wide <- reshape2::dcast(subset(jsum, !is.na(cluster)), formula = cluster ~ cond, value.var = "frac", fill = 0)
  
  jsum.mat <- as.matrix(subset(jsum.wide, select = c(Unenriched, Linneg, StemCell)))
  rownames(jsum.mat) <- jsum.wide$cluster
  return(list(jsum = jsum, jsum.mat = jsum.mat, jsum.wide = jsum.wide))
}



# Check simplex -----------------------------------------------------------

jmark1 <- "H3K4me1"
jmark2 <- "H3K4me3"
jmark3 <- "H3K27me3"
jmark4 <- "H3K9me3"
inf1 <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark1, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
load(inf1, v=T)

dat.umap.glm.fillNAs1 <- dat.umap.glm.fillNAs
dat.umap.glm1 <- dat.umap.glm
dat.umap.lda1 <- dat.umap.lda
mm.celltype.lst1 <- mm.celltype.lst

inf2 <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark2, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
load(inf2, v=T)

dat.umap.glm.fillNAs2 <- dat.umap.glm.fillNAs
dat.umap.glm2 <- dat.umap.glm
dat.umap.lda2 <- dat.umap.lda
mm.celltype.lst2 <- mm.celltype.lst

inf3 <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark3, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
load(inf3, v=T)

dat.umap.glm.fillNAs3 <- dat.umap.glm.fillNAs
dat.umap.glm3 <- dat.umap.glm
dat.umap.lda3 <- dat.umap.lda
mm.celltype.lst3 <- mm.celltype.lst

inf4 <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark4, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
load(inf4, v=T)

dat.umap.glm.fillNAs4 <- dat.umap.glm.fillNAs
dat.umap.glm4 <- dat.umap.glm
dat.umap.lda4 <- dat.umap.lda
mm.celltype.lst4 <- mm.celltype.lst


dat.umap.lda1 <- dat.umap.lda1 %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"))

dat.umap.lda2 <- dat.umap.lda2 %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"))

dat.umap.lda3 <- dat.umap.lda3 %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"))

dat.umap.lda4 <- dat.umap.lda4 %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "_"))

dat.umap.lda1 <- dat.umap.lda1 %>%
  rowwise() %>%
  mutate(cond = GetCondFromSamp(cell, mark = jmark1))

dat.umap.lda2 <- dat.umap.lda2 %>%
  rowwise() %>%
  mutate(cond = GetCondFromSamp(cell, mark = jmark2))

dat.umap.lda3 <- dat.umap.lda3 %>%
  rowwise() %>%
  mutate(cond = GetCondFromSamp(cell, mark = jmark3))

dat.umap.lda4 <- dat.umap.lda4 %>%
  rowwise() %>%
  mutate(cond = GetCondFromSamp(cell, mark = jmark4))
  
m1 <- ggplot(dat.umap.lda1, aes(x = cluster, group = cond, fill = cond)) + geom_bar(position = "dodge") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
m2 <- ggplot(dat.umap.lda2, aes(x = cluster, group = cond, fill = cond)) + geom_bar(position = "dodge") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

m3 <- ggplot(dat.umap.lda3, aes(x = cluster, group = cond, fill = cond)) + geom_bar(position = "dodge") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

m4 <- ggplot(dat.umap.lda4, aes(x = cluster, group = cond, fill = cond)) + geom_bar(position = "dodge") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

multiplot(m1, m2, m3, m4, cols = 2)
  
# Plot simplex  -----------------------------------------------------------

# need to have same number of cells for each condition

jseed <- 0
jsum.out1 <- CalculateCelltypeFractions(dat.umap.lda1, jfrac = 0.8, jseed = jseed)
jsum.out2 <- CalculateCelltypeFractions(dat.umap.lda2, jfrac = 0.8, jseed = jseed)
jsum.out3 <- CalculateCelltypeFractions(dat.umap.lda3, jfrac = 0.8, jseed = jseed)
jsum.out4 <- CalculateCelltypeFractions(dat.umap.lda4, jfrac = 0.8, jseed = jseed)

jsum.mat1 <- jsum.out1$jsum.mat
jsum.mat2 <- jsum.out2$jsum.mat
jsum.mat3 <- jsum.out3$jsum.mat
jsum.mat4 <- jsum.out4$jsum.mat

jsum.mat.merge <- do.call(rbind, list(jsum.mat1, jsum.mat2, jsum.mat3, jsum.mat4))

pca.out <- prcomp(jsum.mat.merge, center = TRUE, scale. = TRUE)
dat.pca <- data.frame(cluster = rownames(pca.out$x), pc1 = pca.out$x[, 1], pc2 = pca.out$x[, 2], stringsAsFactors = FALSE)

head(jsum.out1$jsum.wide)

jsum.out1$jsum.wide$mark <- jmark1
jsum.out2$jsum.wide$mark <- jmark2
jsum.out3$jsum.wide$mark <- jmark3
jsum.out4$jsum.wide$mark <- jmark4


jsum.wide.merge <- do.call(rbind, list(subset(jsum.out1$jsum.wide, select = c(cluster, mark)), 
                                       subset(jsum.out2$jsum.wide, select = c(cluster, mark)),
                                       subset(jsum.out3$jsum.wide, select = c(cluster, mark)),
                                       subset(jsum.out4$jsum.wide, select = c(cluster, mark))))

dat.pca <- left_join(dat.pca, jsum.wide.merge)

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC78A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#eb9d01", "#7fbedf", "#009E73")
ggplot(dat.pca, aes(x = pc1, y = pc2, label = cluster, color = mark)) + geom_point() + 
  theme_bw(24) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_text_repel() + scale_color_manual(values = cbPalette) 



