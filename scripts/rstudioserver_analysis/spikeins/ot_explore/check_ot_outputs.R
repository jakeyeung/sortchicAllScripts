# Jake Yeung
# Date of Creation: 2020-11-24
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/ot_explore/check_ot_outputs.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Check output of ot ------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmark1 <- "H3K4me1"
jmark2 <- "H3K4me3"
inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/glmpca_outputs_factors_for_integration/scotoutputs/scot_output.", jmark1, ".", jmark2, "_coupling.again.txt"))

mat <- as.matrix(fread(inf, header = FALSE, sep = "\t"))



# Load annotations --------------------------------------------------------

inf.jnames1 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_outputs_factors_for_integration/glmpca_factors_", jmark1, ".niter_500.binskeep1000.rownames.txt")
inf.jnames2 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_outputs_factors_for_integration/glmpca_factors_", jmark2, ".niter_500.binskeep1000.rownames.txt")

jnames1 <- fread(inf.jnames1, header = FALSE)
jnames2 <- fread(inf.jnames2, header = FALSE)

rownames(mat) <- unlist(jnames1)
colnames(mat) <- unlist(jnames2)



# Get cell annotations ----------------------------------------------------

jmarks <- c(jmark1, jmark2)
names(jmarks) <- jmarks
indir.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins"
jdate <- "2020-11-18"
dat.metas.init <- lapply(jmarks, function(jmark){
  fname <- paste0("cell_cluster_table_with_spikeins.", jmark, ".", jdate, ".dupfilt.txt")
  inf.meta <- file.path(indir.meta, fname)
  dat.tmp <- fread(inf.meta)
  dat.tmp <- subset(dat.tmp, select = -c(umap1, umap2))
})

library(hash)
library(JFuncs)
library(scchicFuncs)
cell2clusterhash2 <- hash::hash(dat.metas.init[[jmark2]]$cell, dat.metas.init[[jmark2]]$cluster)
cell2clusterhash1 <- hash::hash(dat.metas.init[[jmark1]]$cell, dat.metas.init[[jmark1]]$cluster)

mat.renamed <- mat
colnames(mat.renamed) <- sapply(colnames(mat.renamed), function(x) AssignHash(x = x, jhash = cell2clusterhash2, null.fill = x))

# check a granu 
jcells <- subset(dat.metas.init$H3K4me1, cluster == "Bcells" & cell %in% unlist(jnames1))$cell
names(jcells) <- jcells


# mat.renamed[jcells, ]

# check some cells
jtops.lst <- lapply(jcells, function(jcell){
  jtops <- head(sort(mat.renamed[jcell, ], decreasing = TRUE))
})


tophits <- lapply(jtops.lst, function(x) unique(names(x)))


# Check inputs are correct ------------------------------------------------


niter <- "500"
binskeep <- 1000
jsuffix <- paste0("bincutoff_0.binskeep_", binskeep, ".byplate.szname_none.niter_", niter, ".reorder_rownames.dupfilt")
infs.lst <- lapply(jmarks, function(jmark){
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/glmpca.", jmark, ".", jsuffix, ".RData"))
  assertthat::assert_that(file.exists(inf))
  return(inf)
})


glmpca.factors.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- infs.lst[[jmark]]
  # inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/glmpca.", jmark, ".", jsuffix, ".RData"))
  load(inf, v=T)
  return(glm.out$factors)
})


# Project with matrix -----------------------------------------------------

mat.colnorm <- sweep(x = mat, MARGIN = 2, STATS = colSums(mat), FUN = "/")
mat.rownorm <- sweep(x = mat, MARGIN = 1, STATS = rowSums(mat), FUN = "/")

jproj.colnorm <- t(mat.colnorm) %*% as.matrix(glmpca.factors.lst[[jmark1]])

jproj.rownorm <- mat.rownorm %*% as.matrix(glmpca.factors.lst[[jmark2]])

mat.merge <- rbind(as.matrix(glmpca.factors.lst[[jmark1]]), jproj.rownorm)

library(hash)
library(igraph)
library(umap)
jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

dat.umap.proj <- DoUmapAndLouvain(mat.merge, jsettings)

dat.umap.proj$cluster <- sapply(dat.umap.proj$cell, function(x) AssignHash(x = x, jhash = cell2clusterhash1, null.fill = x))
dat.umap.proj$batch <- rep(c(1, 2), each = nrow(dat.umap.proj) / 2)

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
ggplot(dat.umap.proj, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  facet_wrap(~batch) + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())







