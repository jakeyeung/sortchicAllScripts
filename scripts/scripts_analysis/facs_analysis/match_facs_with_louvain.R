# Jake Yeung
# Date of Creation: 2019-03-18
# File: ~/projects/scchic/scripts/scripts_analysis/facs_analysis/match_facs_with_louvain.R
# Match louvain clusters with facs

library(dplyr)
library(ggplot2)
library(data.table)
library(ggrepel)

# Load the data first -----------------------------------------------------

load("/Users/yeung/data/scchic/robjs/primetime_objs/four_marks_lda_output.RData", v=T)


# Show FACS data  ---------------------------------------------------------

jmark <- "H3K9me3"
# jmark <- "H3K4me1"
# jmark <- "H3K4me3"
# jmark <- "H3K27me3"

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

m <- ggplot(loadings.long, aes(x = PC1, y = PC2, label = facs.feature)) + 
  geom_point() + geom_text_repel() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + ggtitle(jmark)
print(m)

# Take interesting cells plot them in the UMAP  ---------------------------

# tag some cells
cells <- (samps.long %>% arrange(desc(abs(PC1))) %>% mutate(rnk = seq(length(PC1))) %>% filter(rnk <= 20))$samp

m2 <- ggplot(samps.long %>% mutate(samp = ifelse(samp %in% cells, cells, NA)), aes(x = PC1, y = PC2, label = samp)) + 
  geom_point() + geom_text_repel() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + ggtitle(jmark)
print(m2)

m3 <- ggplot(dat.umap.long.lst[[jmark]] %>% rowwise() %>% mutate(cell.lab = ifelse(cell %in% cells, cells, NA)), aes(x = umap1, y = umap2, color = cell.lab)) + 
  geom_point(alpha = 0.5, size = 3) + 
  # geom_text() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(jmark)
print(m3)


