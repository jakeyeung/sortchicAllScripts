# Jake Yeung
# Date of Creation: 2022-07-22
# File: ~/projects/scchic/scripts/revision_scripts/revisions2/2-CCA_full_data_repressive_k9me3_cluster_specific_bins_keeptop_bymark.R
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



library(topicmodels)
library(scchicFuncs)


suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option

parser$add_argument('-refmark', metavar='eg k4me1 k4me3 k27me3',
                    help='refmark')
parser$add_argument('-mark', metavar='eg k4me1 k4me3 k27me3',
                    help='mark')
parser$add_argument('-repressfactor', metavar="adj factor", default = 1, type="integer",
                    help='-1 or 1 to match rep mat to act mat')
parser$add_argument('-outdir', metavar='OUTDIR',
                    help='output dir')
parser$add_argument("-v", "--verbose", action="store_true", default=TRUE,
                    help="Print extra output [default]")


jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 30
jsettings[["min_dist"]] <- 0.1
jsettings[["random_state"]] <- 123

args <- parser$parse_args()
jmarkother <- args$mark
jmarkref <- args$refmark
outdir <- args$outdir

jfactor <- args$repressfactor

# jmarkother <- "k4me1"
# jmarkref <- "k9me3"

jmarksall <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarksall) <- jmarksall
jmarksoldall <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksoldall) <- jmarksall

jmarks <- jmarksall[c(jmarkref, jmarkother)]
jmarksold <- jmarksoldall[c(jmarkref, jmarkother)]

jmarkothers <- jmarksall[jmarkother]

keeptop <- 500
# outdir <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/with_k9me3_cluster_specific_bins_keeptop_", keeptop, "_bymark")
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

inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umaps_pcas_with_batch_corrections/umap_metadata_primetime.k4me1.2022-04-21.txt")
dat.meta.act <- fread(inf.meta) %>%
  left_join(., dat.colors.fixed) %>%
  rowwise() %>%
  mutate(colcode = colcodenew)
dat.meta.colors <- subset(dat.meta.act, select = c(ctype.from.LL, colcode))
ctype2col <- hash::hash(dat.meta.colors$ctype.from.LL, dat.meta.colors$colcode)



# Load impute -------------------------------------------------------------



# lda.main <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_TSS_k9_cluster_specific_bins_keeptop_", keeptop)
lda.main <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_TSS_k9_cluster_specific_bins_keeptop_", keeptop, "_with_TSS")

bad.cells <- readr::read_lines("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/CCA_outputs/k4me1_bad_cells.txt")

dat.impute.lst <- lapply(jmarksold, function(jmark){
  # lda.dir <- file.path(lda.main, paste0("lda_outputs.count_mat_allmerged_for_LDA_k9_cluster_specific_bins_keeptop_", keeptop, ".", jmark, ".2022-07-22"))
  # inf.impute <- file.path(lda.dir, paste0("ldaOut.count_mat_allmerged_for_LDA_k9_cluster_specific_bins_keeptop_", keeptop, ".", jmark, ".2022-07-22.Robj"))
  lda.dir <- file.path(lda.main, paste0("lda_outputs.count_mat_allmerged_for_LDA_k9_cluster_specific_bins_keeptop_", keeptop, "_with_TSS.", jmark, ".2022-07-27"))
  inf.impute <- file.path(lda.dir, paste0("ldaOut.count_mat_allmerged_for_LDA_k9_cluster_specific_bins_keeptop_", keeptop, "_with_TSS.", jmark, ".2022-07-27.Robj"))
  load(inf.impute, v=T) # out.lda
  tm.result <- posterior(out.lda) %>%
    AddTopicToTmResult()
  mat <- log2(t(tm.result$topics %*% tm.result$terms))
  print("remove bad cells")
  print(dim(mat))
  cells.keep <- !colnames(mat) %in% bad.cells
  mat <- mat[, cells.keep]
  print("remove bad cells.. done")
  print(dim(mat))
  mat <- jfactor * mat
  # if (jmark == "k27me3" | jmark == "k9me3"){
  # } else {
  #   mat <- 1 * mat
  # }
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


# Run Liger ---------------------------------------------------------------


# outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/liger_outputs"
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



