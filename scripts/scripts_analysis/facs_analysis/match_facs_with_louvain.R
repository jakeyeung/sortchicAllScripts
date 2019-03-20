# Jake Yeung
# Date of Creation: 2019-03-18
# File: ~/projects/scchic/scripts/scripts_analysis/facs_analysis/match_facs_with_louvain.R
# Match louvain clusters with facs

rm(list=ls())

library(dplyr)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ggrastr)

# Load the data first -----------------------------------------------------

load("/Users/yeung/data/scchic/robjs/primetime_objs/four_marks_lda_output.RData", v=T)

source("scripts/Rfunctions/PlotFunctions.R")

ShortenCellName <- function(x){
  # "BM_H3K9me3_m1_rep1_cell1" -> "m1_rep1_cell1"
  if (is.na(x)){
    return(x)
  }
  return(paste(strsplit(x, "_")[[1]][3:5], collapse = "_"))
}

# Show FACS data  ---------------------------------------------------------

# jmark <- "H3K9me3"
# jmark <- "H3K4me1"
# jmark <- "H3K4me3"
# jmark <- "H3K27me3"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")

jtopn <- 15
pdf(paste0("~/Dropbox/scCHiC_figs/FIG4_BM/facs_integration/facs_integration.chooseTop.", jtopn, ".pdf"), useDingbats = FALSE)
for (jmark in jmarks){
  inf <- file.path(paste0('/Users/yeung/data/scchic/facs/', jmark, '_index_2mice_4plates.csv'))
  dat.facs <- read.table(inf)
  # PCA on the dat.facs
  
  # remove bad columns
  cvar <- apply(dat.facs, 2, var)
  cols.remove <- cvar == 0
  dat.facs <- dat.facs[, !cols.remove]
  
  pca.facs <- prcomp(dat.facs, center = TRUE, scale. = TRUE)
  
  loadings.long <- data.frame(pca.facs$rotation, facs.feature = rownames(pca.facs$rotation), stringsAsFactors = FALSE)
  samps.long <- data.frame(pca.facs$x, samp = rownames(pca.facs$x), stringsAsFactors = FALSE)
  
  
  # plot output
  
  m1 <- ggplot(loadings.long, aes(x = PC1, y = PC2, label = facs.feature)) + 
    geom_point(alpha = 0.3) + geom_text_repel() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + ggtitle(paste(jmark, "FACS feature loadings"))
  
  # Take interesting cells plot them in the UMAP  ---------------------------
  
  # tag some cells
  cells <- (samps.long %>% arrange(desc(abs(PC1))) %>% mutate(rnk = seq(length(PC1))) %>% filter(rnk <= jtopn))$samp
  
  # cells <- (samps.long %>% arrange(desc(PC3)) %>% mutate(rnk = seq(length(samp))) %>% filter(rnk <= 20))$samp
  
  # cells <- (samps.long %>% arrange(PC5) %>% mutate(rnk = seq(length(samp))) %>% filter(rnk <= 20))$samp
  # cells <- (samps.long %>% arrange(desc(PC3)) %>% mutate(rnk = seq(length(samp))) %>% filter(rnk <= 20))$samp
  # cells <- (samps.long %>% arrange(desc(PC3)) %>% mutate(rnk = seq(length(samp))) %>% filter(rnk <= 20))$samp
  
  m2 <- ggplot(samps.long %>% mutate(samp = ifelse(samp %in% cells, cells, NA)) %>% rowwise() %>% mutate(samp = ShortenCellName(samp), samp.bin = ifelse(is.na(samp), "NotSelected", "Selected")), aes(x = PC1, y = PC2, label = samp, color = samp.bin)) + 
    geom_point(alpha = 0.3) + geom_text_repel(color = "black") + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + ggtitle(paste(jmark, "FACS sample loadings. Top:", jtopn)) 
  
  multiplot(m2, m1, cols = 2)
  
  print(m2)
  print(m1)
  
  # m2 <- ggplot(samps.long %>% mutate(samp = ifelse(samp %in% cells, cells, NA)), aes(x = PC1, y = PC3, label = samp)) + 
  #   geom_point() + geom_text_repel() + 
  #   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  #   geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + ggtitle(jmark)
  # print(m2)
  
  jtmp <- RankOrder(dat.umap.long.lst[[jmark]] %>% rowwise() %>% mutate(cell.lab = ifelse(cell %in% cells, "Selected", "NotSelected")), cname = "cell.lab", out.cname = "cell.lab.rank")
  m3 <- ggplot(jtmp, aes(x = umap1, y = umap2, color = cell.lab, order = cell.lab.rank)) + 
    geom_point(alpha = 0.7, size = 3) + 
    # geom_text() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(paste(jmark, "Top:", jtopn))
  print(m3)
}
dev.off()