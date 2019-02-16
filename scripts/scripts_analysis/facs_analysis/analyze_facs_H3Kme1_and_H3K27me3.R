# Jake Yeung
# Date of Creation: 2019-02-16
# File: ~/projects/scchic/scripts/scripts_analysis/facs_analysis/analyze_facs_H3Kme1_and_H3K27me3.R
# Explore both H3K4me1 and H3K27me3

# rm(list=ls())

library(data.table)
library(dplyr)
library(ggplot2)

library(umap)


# Load data ---------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K27me3")

infs.lda <- list("~/Dropbox/scCHiC_figs/FIG4_BM/tables/cluster_tables_for_chloe/cluster_table_H3K4me1.rds",
                 "~/Dropbox/scCHiC_figs/FIG4_BM/tables/cluster_tables_for_chloe/cluster_table_H3K27me3.rds")
names(infs.lda) <- jmarks
dat.inlst <- lapply(infs.lda, function(inf) readRDS(inf))
# label louvain with the mark

for (jmark in jmarks){
  # print(jmark)
  dat.inlst[[jmark]]$mark <- jmark
}
  
dat.inlst <- lapply(dat.inlst, function(dat){
  dat$louvain.mark <- paste(dat$mark, dat$louvain, sep = "-")
  return(dat)
})

dat.long <- bind_rows(dat.inlst)

infs <- list("~/data/scchic/facs/k4me1_index_2mice_4plates.csv", 
             "~/data/scchic/facs/k27me3_index_2mice_4plates.csv")
facs.lst <- lapply(infs, function(inf) fread(inf))


# Batch effects? ----------------------------------------------------------

facs.merge <- bind_rows(facs.lst)

# scale
rnames <- facs.merge$V1
rownames(facs.merge) <- facs.merge$V1
facs.merge$V1 <- NULL

keep.cols.i <- apply(facs.merge, 2, function(jcol) var(jcol)) > 0

facs.merge <- as.matrix(facs.merge)[, keep.cols.i]

facs.merge <- scale(facs.merge, center = TRUE, scale = TRUE)

rownames(facs.merge) <- rnames

custom.settings <- umap::umap.defaults
custom.settings$n_neighbors <- 40
custom.settings$metric <- "euclidean"
custom.settings$min_dist <- 0.2
custom.settings$random_state <- 0

umap.out <- umap::umap(facs.merge, config = custom.settings)

colvec <- as.numeric(as.factor(sapply(rnames, function(x) strsplit(x, "_")[[1]][[2]])))

plot(umap.out$layout[, 1], umap.out$layout[, 2], pch = 20, col = colvec)

# color by louvain in H3K4me1
louvain.hash <- hash(dat.long$cell, dat.long$louvain.mark)

umap.out.long <- data.frame(umap1 = umap.out$layout[, 1],
                            umap2 = umap.out$layout[, 2],
                            cell = rownames(umap.out$layout), stringsAsFactors = FALSE)
umap.out.long$louvain.mark <- sapply(rownames(umap.out$layout), function(x) louvain.hash[[x]])
umap.out.long$mark <- sapply(umap.out.long$cell, function(x) strsplit(x, "_")[[1]][[2]])

umap.out.long <- subset(umap.out.long, !is.null(louvain.mark))

umap.out.long <- umap.out.long %>%
  mutate(col = as.character(ifelse(mark == "H3K4me1", louvain.mark, "zH3K27me3")))
  # mutate(col = as.character(ifelse(louvain.mark %in% c("H3K4me1-5", "H3K4me1-8"), louvain.mark, "other")))

# some missing values? redo this later
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#FFD700")

m1 <- ggplot(umap.out.long, aes(x = umap1, y = umap2, color = as.factor(col))) + geom_point(size=0.5)  + facet_wrap(~mark) + 
  scale_colour_manual(values=cbPalette)+ theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Recall what are the louvain in regulatory space -------------------------

m2 <- ggplot(dat.long %>% filter(mark == "H3K4me1" & louvain != 4), aes(x = umap1, y = umap2, color = as.factor(louvain.mark))) + geom_point(size=0.5) + facet_wrap(~mark) + 
  scale_colour_manual(values=cbPalette) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# print(m2)

multiplot(m1, m2, col = 2)


