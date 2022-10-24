# Jake Yeung
# Date of Creation: 2022-07-27
# File: ~/projects/scchic/scripts/revision_scripts/revisions2/2-Liger_full_data_repressive.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
# library(Seurat)
library(rliger)

library(hash)
library(igraph)
library(umap)
library(scchicFuncs)
jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 30
jsettings[["min_dist"]] <- 0.1
jsettings[["random_state"]] <- 123

jmarks <- c("k4me1", "k4me3", "k27me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarksold) <- jmarks

jmarkref <- jmarks[[2]]
jmarkothers <- jmarks[jmarks != jmarkref]

# outdir <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/liger_outputs/refmark_", jmarkref)
# dir.create(outdir)

# Load new colors ---------------------------------------------------------

inf.colors.fixed <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/ctypes_on_umap_batch_corrected_colors_fixed/dat_colors_DC_monocyte_fixed.2022-05-17.txt"
dat.colors.fixed <- fread(inf.colors.fixed)


# Load meta ----------------------------------------------------------------

dat.meta.lst <- lapply(jmarks, function(jmark){
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umaps_pcas_with_batch_corrections/umap_metadata_primetime.", jmark, ".2022-04-21.txt")
  dat.meta <- fread(inf.meta) %>%
    left_join(., dat.colors.fixed) %>%
    rowwise() %>%
    mutate(colcode = colcodenew)
  # replace colcode with colcodenew
})

dat.meta.colors <- subset(dat.meta.lst$k4me1, select = c(ctype.from.LL, colcode))
ctype2col <- hash::hash(dat.meta.colors$ctype.from.LL, dat.meta.colors$colcode)




# Get marker genes --------------------------------------------------------

dat.markers <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/genesets/dat_genes_sub_join_from_orig.rds")
print(unique(dat.markers$jset))

genes.keep <- unique(dat.markers$gene)


# Load impute -------------------------------------------------------------


dat.impute.lst <- lapply(jmarks, function(jmark){
  inf.impute <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/primetime_objects/dat_impute_tss_", jmark, ".2022-04-24.rds")
  mat <- readRDS(inf.impute)
  if (jmark == "k27me3"){
    mat <- -1 * mat
  } else {
    mat <- 1 * mat
  }
  return(mat)
  # mat[genes.keep, ]
})

dat.impute.markers.lst <- lapply(dat.impute.lst, function(jdat){
  rows.keep <- rownames(jdat) %in% genes.keep
  jdat[rows.keep, ]
})

common.markers <- Reduce(intersect, lapply(dat.impute.markers.lst, function(jdat) rownames(jdat)) )

dat.meta.lst <- lapply(jmarks, function(jmark){
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umaps_pcas_with_batch_corrections/umap_metadata_primetime.", jmark, ".2022-04-21.txt")
  dat.meta <- fread(inf.meta) %>%
    left_join(., dat.colors.fixed) %>%
    rowwise() %>%
    mutate(colcode = colcodenew) %>%
    as.data.frame()
  rownames(dat.meta) <- dat.meta$cell
  return(dat.meta)
  # replace colcode with colcodenew
})


dat.umap.annot.unif <- lapply(jmarks, function(jmark){
  jdat <- dat.meta.lst[[jmark]]
  jdat.select <- subset(jdat, select = c(cell, ctype.from.LL, colcode))
  jdat.select$mark <- jmark
  return(jdat.select)
}) %>%
  bind_rows()


# Set up liger objects ----------------------------------------------------


# seurat.lst <- lapply(jmarks, function(jmark){
#   print(jmark)
#   chic <- CreateSeuratObject(counts = exp(dat.impute.markers.lst[[jmark]]), meta.data = dat.meta.lst[[jmark]], assay = jmark)
#   chic <- NormalizeData(chic)
#   chic <- ScaleData(chic, do.scale = TRUE)
#   chic <- FindVariableFeatures(chic, nfeatures = 500, verbose = TRUE)
#   chic@meta.data$mark <- jmark
#   return(chic)
# })

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/liger_outputs"
marks.combined.lst <- lapply(jmarkothers, function(jmark){
  print(jmark)
  raw.data <- list(ref = exp(dat.impute.markers.lst[[jmarkref]]), target =  exp(dat.impute.markers.lst[[jmark]]))
  common.genes <- intersect(rownames(raw.data[[1]]), rownames(raw.data[[2]]))
  names(raw.data) <- c(jmarkref, jmark)
  
  ligerobj <- rliger::createLiger(raw.data = raw.data)
  ligerobj <- rliger::normalize(ligerobj)
  # ligerobj <- rliger::selectGenes(ligerobj, var.thresh = -100, alpha.thresh = 1000, tol = 0.001, combine = "union", num.genes = 500, do.plot = TRUE)
  ligerobj@var.genes <- common.genes
  ligerobj <- rliger::scaleNotCenter(ligerobj, remove.missing = FALSE)
  ligerobj <- rliger::optimizeALS(ligerobj, k = 20)
  
  saveRDS(ligerobj, file = file.path(outdir, paste0("liger_", jmarkref, "_and_", jmark, "_optimizeALS.", Sys.Date(), ".rds")))
  ligerobj <- rliger::quantile_norm(ligerobj)
  ligerobj <- rliger::louvainCluster(ligerobj, resolution = 0.2)
  ligerobj <- rliger::runUMAP(ligerobj, distance = 'cosine', n_neighbors = 30, min_dist = 0.3)
  saveRDS(ligerobj, file = file.path(outdir, paste0("liger_", jmarkref, "_and_", jmark, "_UMAP.", Sys.Date(), ".rds")))
  return(ligerobj)
})
saveRDS(marks.combined.lst, file = file.path(outdir, paste0("liger_all_marks_combined.", Sys.Date(), ".rds")))