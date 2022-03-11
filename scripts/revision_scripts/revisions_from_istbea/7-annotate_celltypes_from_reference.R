# Jake Yeung
# Date of Creation: 2022-02-01
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/7-annotate_celltypes_from_reference.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(topicmodels)
library(AnnotateCelltypes)
library(parallel)

jmarks <- c("k4me1", "k4me3", "k27me3"); names(jmarks) <- jmarks

jmark <- jmarks[[2]]

LL.all.lst <- mclapply(jmarks, function(jmark){
  
  
  # Load metadata -----------------------------------------------------------
  
  indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/varfilt"
  inf.meta <- file.path(indir.meta, paste0("metadata_celltyping_", jmark, ".allbins.2022-02-01.txt"))
  dat.meta <- fread(inf.meta)
  
  print(unique(sort(dat.meta$ctype)))
  
  # reclassify celltypes for AllCells, HSPCs, IL7RLinNeg, LinNeg, LSK
  
  
  # Load rawcounts ----------------------------------------------------------
  
  indir.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_allbins.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28")
  inf.ldaout <- file.path(indir.ldaout, paste0("ldaOut.count_mat_var_filt_allbins.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28.Robj"))
  load(inf.ldaout, v=T)
  
  
  
  # Annotate celltypes ------------------------------------------------------
  
  ctypes.to.annotate <- c("Basophils", "Bcells", "CMP", "DCs", "Eryths", "GMP", "Granulocytes", "HSCs", "LT", "MEP", "Monocytes", "MPPs", "NKs", "pDCs", "ST", "Tcells")
  names(ctypes.to.annotate) <- ctypes.to.annotate
  
  dat.meta.ctype.filt <- subset(dat.meta, ctype %in% ctypes.to.annotate)
  dat.meta.testset <- subset(dat.meta, !ctype %in% ctypes.to.annotate)
  cluster.ids <- lapply(split(dat.meta.ctype.filt, f = dat.meta.ctype.filt$ctype), function(x) x$cell)
  
  
  count.mat.pbulk <- SumAcrossClusters(count.mat = count.mat, cnames.keep.lst = cluster.ids)
  count.mat.pbulk <- do.call(cbind, count.mat.pbulk)
  
  
  #' zero counts are problematic when calculating likelihoods, add pseudocount to compensate
  mats.keep.pseudobulk.pcount <- count.mat.pbulk + 1
  mats.keep.pseudobulk.pcount.frac <- sweep(mats.keep.pseudobulk.pcount, MARGIN = 2, STATS = colSums(mats.keep.pseudobulk.pcount), FUN = "/")
  
  
  # Estimate from multionomials  --------------------------------------------
  
  
  # (cellname <- dat.meta.testset$cell[[1]])
  # jcell <- count.mat[, cellname]
  # 
  # #' note the input here is unnormalized counts. No need for log-transform or scaling rows or columns.
  # #' should be robust to downsampling too.
  # LL <- AnnotateCelltypes::CalculateMultinomLL(jcell, mats.keep.pseudobulk.pcount.frac)
  
  
  # run on all
  LL.all <- apply(count.mat, 2, function(jcell){
    LL <- AnnotateCelltypes::CalculateMultinomLL(jcell, mats.keep.pseudobulk.pcount.frac)
  })
  rownames(LL.all) <- colnames(mats.keep.pseudobulk.pcount.frac)
  
  outrds <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping/LLmat_allbins", ".", jmark, ".", Sys.Date(), ".rds")
  saveRDS(LL.all, file = outrds)
  
  
}, mc.cores = length(jmarks))




# saveoutputs
# LL.mat <- t(LL.all)

# 
# 
# # Check best hits ---------------------------------------------------------
# 
# ctypes.indx <- colnames(LL.mat)
# 
# LL.best <- apply(LL.mat, 1, function(jrow){
#   ctypes.indx[which.max(jrow)]
# })
# 
# dat.LL.best <- data.frame(cell = names(LL.best), ctype.from.LL = LL.best, stringsAsFactors = FALSE)
# 
# dat.meta.reannotate <- dat.meta %>%
#   left_join(., dat.LL.best)
# 
# dat.meta.reannotate.sum <- dat.meta.reannotate %>%
#   group_by(ctype, ctype.from.LL) %>%
#   summarise(ncells = length(cell)) %>%
#   group_by(ctype) %>%
#   mutate(nfracs = ncells / sum(ncells))
# 
# 
# print(sort(unique(dat.meta.reannotate.sum$ctype)))
# print(sort(unique(dat.meta.reannotate.sum$ctype.from.LL)))
# 
# dat.meta.reannotate.sum2 <- dat.meta.reannotate %>%
#   group_by(ctype, ctype.from.LL) %>%
#   summarise(ncells = length(cell)) %>%
#   group_by(ctype.from.LL) %>%
#   mutate(nfracs = ncells / sum(ncells))
# 
# 
# 
# # Show cleaned UMAP  ------------------------------------------------------
# 
# library(hash)
# library(igraph)
# library(umap)
# 
# jsettings <- umap.defaults
# jsettings[["n_neighbors"]] <- 30
# jsettings[["min_dist"]] <- 0.1
# jsettings[["random_state"]] <- 123
# 
# pca.LL <- prcomp(LL.all, center = TRUE, scale. = TRUE)
# jsettings <- umap.defaults
# jsettings$n_neighbors <- 30
# jsettings$min_dist <- 0.1
# jsettings$random_state <- 123
# 
# 
# umap.out <- umap(pca.LL$rotation, config = jsettings)
# dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
# dat.umap.long.annot <- left_join(dat.umap.long, subset(dat.meta.reannotate, select = c(-umap1, -umap2)))
# 
# cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#d9d331")
# 
# ggplot(dat.umap.long.annot, aes(x = umap1, y = umap2, color = ctype.from.LL)) +
#   geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   scale_color_manual(values = cbPalette) +
#   ggtitle("Cell type predictions from cluster-free model", jmark)
# 
# 
# ggplot(dat.umap.long.annot, aes(x = umap1, y = umap2, color = ctype.from.LL)) +
#   geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   facet_wrap(~ctype.from.LL) + 
#   scale_color_manual(values = cbPalette) +
#   ggtitle("Cell type predictions from cluster-free model", jmark)
# 
# 
# 
