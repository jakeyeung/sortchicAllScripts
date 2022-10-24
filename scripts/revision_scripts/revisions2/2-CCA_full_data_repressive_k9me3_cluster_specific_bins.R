# Jake Yeung
# Date of Creation: 2022-07-22
# File: ~/projects/scchic/scripts/revision_scripts/revisions2/2-CCA_full_data_repressive_k9me3_cluster_specific_bins.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(Seurat)

library(hash)
library(igraph)
library(umap)
library(scchicFuncs)
jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 30
jsettings[["min_dist"]] <- 0.1
jsettings[["random_state"]] <- 123


jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksold) <- jmarks

jmarkref <- "k9me3"
jmarkothers <- jmarks[jmarks != jmarkref]

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/with_k9me3_cluster_specific_bins"
dir.create(outdir)

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



# # Get marker genes --------------------------------------------------------
# 
# dat.markers <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/genesets/dat_genes_sub_join_from_orig.rds")
# print(unique(dat.markers$jset))
# genes.keep <- unique(dat.markers$gene)



# Load impute -------------------------------------------------------------


lda.main <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_TSS_k9_cluster_specific_bins"

library(topicmodels)
library(scchicFuncs)
dat.impute.lst <- lapply(jmarksold, function(jmark){
  lda.dir <- file.path(lda.main, paste0("lda_outputs.count_mat_allmerged_for_LDA_k9_cluster_specific_bins.", jmark, ".2022-07-21"))
  inf.impute <- file.path(lda.dir, paste0("ldaOut.count_mat_allmerged_for_LDA_k9_cluster_specific_bins.", jmark, ".2022-07-21.Robj"))
  load(inf.impute, v=T) # out.lda
  tm.result <- posterior(out.lda) %>%
    AddTopicToTmResult()
  mat <- log2(t(tm.result$topics %*% tm.result$terms))
  if (jmark == "k27me3" | jmark == "k9me3"){
    mat <- -1 * mat
  } else {
    mat <- 1 * mat
  }
  return(mat)
  # mat[genes.keep, ]
})

common.markers <- Reduce(intersect, lapply(dat.impute.lst, function(jdat) rownames(jdat)) )

dat.impute.markers.lst <- lapply(dat.impute.lst, function(jdat){
  jdat[common.markers, ]
})


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

# Use CCA  ----------------------------------------------------------------


seurat.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  chic <- CreateSeuratObject(counts = exp(dat.impute.markers.lst[[jmark]]), meta.data = dat.meta.lst[[jmark]], assay = jmark)
  chic <- NormalizeData(chic)
  chic <- ScaleData(chic, do.scale = TRUE)
  chic <- FindVariableFeatures(chic, nfeatures = 500, verbose = TRUE)
  chic@meta.data$mark <- jmark
  return(chic)
})


marks.combined.lst <- lapply(jmarkothers, function(jmark){
  print(jmark)
  marks.combined <- RunCCA(seurat.lst[[jmarkref]], seurat.lst[[jmark]], num.cc = 30, assay1 = jmarkref, assay2 = jmark)
  save(marks.combined, file = file.path(outdir, paste0("CCA_", jmarkref, "_and_", jmark, ".repressive.", Sys.Date(), ".RData")))
  return(marks.combined)
})


# Check downstream  -------------------------------------------------------

# jmarkcheck <- jmarks[[3]]

print("Downstream UMAPs")
for (jmarkcheck in jmarkothers){
  
  print(jmarkcheck)
  
  marks.combined <- marks.combined.lst[[jmarkcheck]]
  
  dat.cca <- marks.combined@reductions$cca@cell.embeddings
  
  # loadings <- 
  
  print(dim(dat.cca))
  
  
  dat.umap.cca <- scchicFuncs::DoUmapAndLouvain(dat.cca, jsettings)
  
  save(dat.umap.cca, file = file.path(outdir, paste0("UMAP_of_CCA_repressive_", jmarkref, "_and_", jmarkcheck, ".", Sys.Date(), ".RData")))
  
  dat.umap.cca.annot <- dat.umap.cca %>%
    left_join(., dat.umap.annot.unif, by = "cell")
  
  cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597", "#ff9f7d", "#eb9d01", "#7fbedf")
  
  dat.meta.colors <- subset(dat.meta.lst$k4me1, select = c(ctype.from.LL, colcode))
  
  pdf(file.path(outdir, paste0("plots_CCA_", jmarkref, "_and_", jmarkcheck, "_repressive.pdf")), useDingbats = FALSE)
  m.all <- ggplot(dat.umap.cca.annot, aes(x = umap1, y = umap2, color = mark)) +
    geom_point(alpha = 0.5) +
    xlab("CCA-UMAP1") + ylab("CCA-UMAP2") +
    scale_color_manual(values = cbPalette) +
    theme_bw() +
    ggtitle(jmarkcheck) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m.all)
  
  m.ctype <- ggplot(dat.umap.cca.annot, aes(x = umap1, y = umap2, color = colcode)) +
    geom_point(alpha = 0.5) +
    ggtitle("Colored by dataset") +
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode, guide = "legend") +
    ggtitle(jmarkcheck) + 
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m.ctype)
  dev.off()
  
  
  
}

print("Done")

