# Jake Yeung
# Date of Creation: 2022-07-21
# File: ~/projects/scchic/scripts/revision_scripts/revisions2/2-CCA_full_data_repressive_k9me3_split1.R
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


jmarks <- c("k4me1", "k9me3", "k4me3", "k27me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K9me3", "H3K4me3", "H3K27me3"); names(jmarksold) <- jmarks

jmarkref <- jmarks[[1]]
jmarkothers <- jmarks[jmarks != jmarkref]

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/with_k9me3"
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







# Get marker genes --------------------------------------------------------

dat.markers <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/genesets/dat_genes_sub_join_from_orig.rds")
print(unique(dat.markers$jset))

genes.keep <- unique(dat.markers$gene)



# Load impute -------------------------------------------------------------


lda.main <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_TSS_k9_bins"

library(topicmodels)
library(scchicFuncs)
dat.impute.lst <- lapply(jmarksold, function(jmark){
  lda.dir <- file.path(lda.main, paste0("lda_outputs.count_mat_allmerged_for_LDA.", jmark, ".2022-07-20"))
  inf.impute <- file.path(lda.dir, paste0("ldaOut.count_mat_allmerged_for_LDA.", jmark, ".2022-07-20.Robj"))
  # mat <- readRDS(inf.impute)
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


dat.impute.markers.lst <- lapply(dat.impute.lst, function(jdat){
  rows.keep1 <- rownames(jdat) %in% genes.keep
  rows.keep2 <- rownames(jdat)[!grepl(";NM_", rownames(jdat))]
  rows.keep.cat <- unique(c(rows.keep1, rows.keep2))
  rows.keep.all <- rownames(jdat)[rownames(jdat) %in% rows.keep.cat]
  jdat[rows.keep.all, ]
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

for (jmark in jmarkothers){
  print(jmark)
  marks.combined <- RunCCA(seurat.lst[[jmarkref]], seurat.lst[[jmark]], num.cc = 30, assay1 = jmarkref, assay2 = jmark)
  save(marks.combined, file = file.path(outdir, paste0("CCA_", jmarkref, "_and_", jmark, ".repressive1.", Sys.Date(), ".RData")))
}


