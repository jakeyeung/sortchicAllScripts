# Jake Yeung
# Date of Creation: 2019-02-11
# File: ~/projects/scchic/scripts/scripts_analysis/facs_analysis/analyze_facs_data.R
# FACS analysis


library(dplyr)
library(ggplot2)
library(data.table)
library(umap)
library(hash)

# Load data ---------------------------------------------------------------

# dat <- fread("~/data/scchic/facs/k4me1_index_2mice_4plates.csv", sep = "\t")
dat <- read.table("~/data/scchic/facs/k4me1_index_2mice_4plates.csv", sep = "\t", stringsAsFactors = FALSE)

dat.lda <- readRDS("/Users/yeung/Dropbox/scCHiC_figs/FIG4_BM/tables/cluster_tables_for_chloe/cluster_table_H3K4me1.rds")


# Umap the data -----------------------------------------------------------

dat.scaled <- scale(dat, center = TRUE, scale = TRUE)
dat.scaled <- t(scale(t(dat), center = TRUE, scale = TRUE))

dat.umap <- umap(dat.scaled)

plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20)

dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], 
                            umap2 = dat.umap$layout[, 2], 
                            cell = rownames(dat.umap$layout), stringsAsFactors = FALSE)

# Color by louvain  -------------------------------------------------------

jhash <- hash(dat.lda$cell, dat.lda$louvain)

dat.umap.long$louvain <- sapply(dat.umap.long$cell, function(x){
  xout <- jhash[[x]]
  if (is.null(xout)){
    return(NA)
  } else {
    xout
  }
})

jsub <- subset(dat.umap.long, !is.na(louvain))

ggplot(jsub, aes(x = umap1, y = umap2, color = as.factor(louvain))) + 
  geom_point() + 
  scale_color_brewer(palette = "Spectral") + 
  theme_bw() 

ggplot(dat.lda, aes(x = umap1, y = umap2, color = as.factor(louvain))) + geom_point() + 
  theme_bw() + scale_color_brewer(palette = "Spectral")

ggplot(dat.lda %>% mutate(louvain2=ifelse(louvain==1, TRUE, FALSE)), aes(x = umap1, y = umap2, color = as.factor(louvain2))) + 
  geom_point() + 
  theme_bw() + scale_color_brewer(palette = "Spectral")

ggplot(dat.lda %>% mutate(louvain2=ifelse(louvain %in% c("2", "4", "7", "6"), TRUE, FALSE)), aes(x = umap1, y = umap2, color = as.factor(louvain2))) + 
  geom_point() + 
  theme_bw() + scale_color_brewer(palette = "Spectral")

ggplot(dat.lda %>% mutate(louvain2=ifelse(louvain %in% c("3", "8", "5"), TRUE, FALSE)), aes(x = umap1, y = umap2, color = as.factor(louvain2))) + 
  geom_point() + 
  theme_bw() + scale_color_brewer(palette = "Spectral")

# merge some louvains together
keys <- as.character(seq(8))
vals <- c("1", "2,4,6,7", "3,5,8", "2,4,6,7", "3,5,8", "2,4,6,7", "2,4,6,7", "3,5,8")

clstr2 <- hash(keys, vals)

jsub$louvain2 <- sapply(as.character(jsub$louvain), function(x) clstr2[[x]])

ggplot(jsub, aes(x = umap1, y = umap2, color = as.factor(louvain2))) + 
  geom_point() + 
  scale_color_brewer(palette = "RdYlBu") + 
  theme_bw() 

