# Jake Yeung
# Date of Creation: 2022-04-20
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/24-check_batch_effect_downstream_k27me3.R
# 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 60
jsettings[["min_dist"]] <- 0.5
jsettings[["random_state"]] <- 123
jsettings[["spread"]] <- 4



# Load mat ----------------------------------------------------------------

inf.mat <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/mat_wide_k27me3_batch_corrected.2022-04-19.rds"
mat <- readRDS(inf.mat)

# SVD ---------------------------------------------------------------------

ntopics <- 30

# logodds <- mat  # alreay in -Inf Inf
# # remove mean and SVD
# logodds.centered <- t(scale(t(logodds), center = TRUE, scale = FALSE))
# # logodds.centered.check <- sweep(logodds, MARGIN = 1, STATS = rowMeans(logodds), FUN = "-")
# logodds.pca <- prcomp(t(logodds.centered), center = FALSE, scale. = FALSE, rank. = ntopics)

logodds.pca <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/pca_output_mat_adj.2022-04-20.rds")

PoV <- signif(logodds.pca$sdev^2/sum(logodds.pca$sdev^2), digits = 2)
# summary(logodds.pca)

U.init <- logodds.pca$x  # cells by k
V.init <- logodds.pca$rotation  # genes by k, no need to transpose


# saveRDS(logodds.pca, file = paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/pca_output_mat_adj.", Sys.Date(), ".rds"))

# UMAP --------------------------------------------------------------------


indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths2"
inf.meta <- file.path(indir.meta, paste0("metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.k27me3.txt"))
dat.meta <- fread(inf.meta) %>%
  dplyr::select(cell, ctype.from.LL, batch, colcode) %>%
  rowwise() %>%
  mutate(cluster = ctype.from.LL,
         jrep2 = ifelse(batch == "New", "anew", "old"),
         cluster = ifelse(cluster == "HSCs", "aHSCs", cluster))


dat.meta.colors <- subset(dat.meta, select = c(colcode, ctype.from.LL))
dat.meta.colors <- dat.meta.colors[!duplicated(dat.meta.colors), ]

dat.umap.annot <- DoUmapAndLouvain(U.init, jsettings) %>%
  left_join(., dat.meta)

# dat.umap.annot <- left_join(dat.umap, dat.meta.colors)

m <- ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = batch)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~batch) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
print(m)


m <- ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = colcode)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                        guide = "legend") + 
  facet_wrap(~batch) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
print(m)


dat.pca <- data.frame(cell = rownames(U.init), pc1 = U.init[, 1], pc2 =U.init[, 2], stringsAsFactors = FALSE) %>%
  left_join(., dat.meta)

ggplot(dat.pca, aes(x = pc1, y = pc2, color = colcode)) + 
  geom_point() + 
  xlab(paste0("PC1 (", PoV[[1]], ")")) + 
  ylab(paste0("PC2 (", PoV[[2]], ")")) + 
  theme_bw() + 
  scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                        guide = "legend") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  


