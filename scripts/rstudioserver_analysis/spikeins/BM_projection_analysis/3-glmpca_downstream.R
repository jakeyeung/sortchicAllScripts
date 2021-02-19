# Jake Yeung
# Date of Creation: 2020-12-28
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_projection_analysis/3-glmpca_downstream.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
jsettings$spread <- 8
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


# Load metas --------------------------------------------------------------

dat.metas <- lapply(jmarks, function(jmark){
  inf.meta <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/count_tables.BMround2.from_peaks.sitecount_mat.split_old_and_new/count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.metadata.2020-12-27.txt"))
  fread(inf.meta)
})

# Load glmpca and check ---------------------------------------------------

jmark <- "H3K4me1"
# jsuffix <- "new_to_old"
jsuffix <- "old_to_new"
niter <- "500"
# binskeep <- "1000"
binskeep <- "0"
platename <- "jrep"
# platename <- "experi"

infglmpca <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/BM_from_projections/glmpca.", jmark, ".bincutoff_0.binskeep_", binskeep, ".platename_", platename, ".szname_none.niter_", niter, ".reorder_rownames.dupfilt.suffix_", jsuffix, ".RData"))
assertthat::assert_that(file.exists(infglmpca))

load(infglmpca, v=T)


# Check UMAP  -------------------------------------------------------------

dat.umap <- DoUmapAndLouvain(glm.out$factors, jsettings)

dat.umap.annot <- left_join(dat.umap, dat.metas[[jmark]] %>% dplyr::select(c(cell, cluster, batch, jrep)))

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point(size = 0.5, alpha = 0.8) + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme_bw() + 
  facet_wrap(~jrep) + 
  ggtitle(jmark) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point(size = 0.5, alpha = 0.8) + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme_bw() + 
  ggtitle(jmark) + 
  facet_wrap(~jrep) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point(size = 0.5, alpha = 0.8) + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme_bw() + 
  ggtitle(jmark) + 
  facet_wrap(~cluster) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point(size = 0.5, alpha = 0.25) + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme_bw() + 
  ggtitle(jmark) + 
  facet_wrap(~batch) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
