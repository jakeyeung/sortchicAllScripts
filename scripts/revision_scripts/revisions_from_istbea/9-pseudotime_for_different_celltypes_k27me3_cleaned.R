# Jake Yeung
# Date of Creation: 2022-04-16
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/9-pseudotime_for_different_celltypes_k27me3_cleaned.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmark <- "k27me3"


# Load metadata -----------------------------------------------------------

inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths2/metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.", jmark, ".2022-04-16.txt")
dat.meta.LL <- fread(inf.meta)


# Load raw counts ---------------------------------------------------------

# inf.glmpca <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/glmpca_outputs/varfilt/glmpca.", jmark, ".bincutoff_0.binskeep_0.platename_batch.szname_none.niter_1000.RData")
# inf.glmpca <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/glmpca_outputs/clean_again_k4me3_k27me3/glmpca.", jmark, ".bincutoff_0.binskeep_0.platename_batch.szname_none.niter_1000.allcells.RData")
# inf.glmpca <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/glmpca_outputs/k9me3_renamed/glmpca.k9me3.bincutoff_0.binskeep_0.platename_batch.szname_none.niter_1000.allcells.RData")
inf.glmpca <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/glmpca_outputs/clean_again_k27me3/glmpca.", jmark, ".bincutoff_0.binskeep_0.platename_batch.szname_none.niter_1000.new_only.RData")

load(inf.glmpca, v=T)
count.mat <- glm.inits$Y.filt

dat.glmpca <- data.frame(cell = rownames(glm.out$factors), glm.out$factors, stringsAsFactors = FALSE) %>%
  left_join(., dat.meta.LL)


# Create different trajectories -------------------------------------------

hpcs <- c("HSCs", "LT", "ST", "MPPs")
progs <- c(hpcs, c("CMP", "GMP"))
# granu tracj
ctypes.granus <- c(hpcs, "CMP", "GMP", "Basophils", "Granulocytes")
dims.granus <- c("dim2", "dim3")

ctypes.bcells <- c(hpcs, "Bcells")
dims.bcells <- c("dim5", "dim6")

ctypes.eryths <- c(hpcs, "MEP", "Eryths")
dims.eryths <- c("dim1", "dim2")

ctypes.pdcs <- c(hpcs, "pDCs")
dims.pdcs <- c("dim1", "dim2")

ctypes.dcs <- c(hpcs, "CMP", "GMP", "DCs")
dims.dcs <-  c("dim2", "dim3")

ctypes.basos <- c(hpcs, "CMP", "GMP", "Basophils")
dims.basos <- c("dim1", "dim2")

ctypes.nks <- c(hpcs, "NKs")
dims.nks <- c("dim5", "dim6")

ctypes.monos <- c(hpcs, "CMP", "GMP", "Monocytes")
dims.monos <- c("dim2", "dim3")

ctypes.lst <- list(ctypes.granus, ctypes.bcells, ctypes.eryths, ctypes.pdcs, ctypes.dcs, ctypes.basos, ctypes.nks, ctypes.monos)
names(ctypes.lst) <- sapply(ctypes.lst, function(x) x[length(x)])

ctypes.names <- names(ctypes.lst); names(ctypes.names) <- ctypes.names

dims.lst  <- list(dims.granus, dims.bcells, dims.eryths, dims.pdcs, dims.dcs, dims.basos, dims.nks, dims.monos)
names(dims.lst) <- ctypes.names

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/cleaned2"
dir.create(outdir)
outpdf <- file.path(outdir, paste0("trajs_plots.", jmark, ".", Sys.Date(), ".pdf"))
outrds <- file.path(outdir, paste0("trajs_outputs.", jmark, ".", Sys.Date(), ".rds"))
pdf(outpdf, useDingbats = FALSE)
dat.traj.lst <- lapply(ctypes.names, function(jctype){
  
  jtitle <- paste(jctype, "trajectory")
  ctypes.vec <- ctypes.lst[[jctype]]
  jdims <- dims.lst[[jctype]]
  jdim1 <- jdims[[1]]
  jdim2 <- jdims[[2]]
  
  # plot on glmpca
  dat.glmpca.keep <- dat.glmpca %>%
    rowwise() %>%
    mutate(is.ctype = ctype.from.LL %in% ctypes.vec)
  
  if (jdim1 == "dim1"){
    print("Removing outliers ")
    print(jctype)
    dat.glmpca.keep <- dat.glmpca.keep %>%
      rowwise() %>%
      mutate(is.ctype = ifelse(dim1 > 20 & ctype.from.LL %in% progs, FALSE, is.ctype))
  } else if (jdim2 == "dim3"){
    print("Removing outliers 2")
    print(jctype)
    dat.glmpca.keep <- dat.glmpca.keep %>%
      rowwise() %>%
      mutate(is.ctype = ifelse(dim3 < -5 & ctype.from.LL %in% progs, FALSE, is.ctype))
  }
  
  dat.glmpca.keep <- dat.glmpca.keep %>%
    arrange(is.ctype)
  
  m <- ggplot(dat.glmpca.keep, aes_string(x = jdim1, y = jdim2, color = "is.ctype")) + 
    geom_point() + 
    theme_bw() + 
    # facet_wrap(~ctype.from.LL) + 
    ggtitle(jtitle) + 
    # geom_abline(slope = 0.5, intercept = -3) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.glmpca.keep, aes_string(x = jdim1, y = jdim2, color = "is.ctype")) + 
    geom_point() + 
    theme_bw() + 
    facet_wrap(~ctype.from.LL) + 
    ggtitle(jtitle) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  dat.glmpca.keep.filt <- dat.glmpca.keep %>%
    rowwise() %>%
    mutate(ctype.from.LL.adj = ifelse(is.ctype, ctype.from.LL, NA)) %>%
    ungroup() %>%
    mutate(ctype.ordered = factor(ctype.from.LL.adj, levels = ctypes.vec))
  
  dat.glmpca.keep.filt <- dat.glmpca.keep.filt[order(dat.glmpca.keep.filt$ctype.ordered, na.last = FALSE), ]
  
  m <- ggplot(dat.glmpca.keep.filt, aes_string(x = jdim1, y = jdim2, color = "ctype.ordered")) + 
    geom_point() + 
    scale_color_viridis_d(na.value = "grey85") + 
    theme_bw() + 
    ggtitle(jtitle) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  # get pseudotime by PCA
  dat.filt <- subset(dat.glmpca.keep.filt, !is.na(ctype.ordered))
  mat.filt <- as.matrix(dat.filt[, jdims]); rownames(mat.filt) <- dat.filt$cell
  pca.filt <- prcomp(mat.filt, center = TRUE, scale. = TRUE)
  
  dat.ptime.filt <- data.frame(cell = rownames(pca.filt$x), ptime = pca.filt$x[, 1], stringsAsFactors = FALSE) %>%
    ungroup() %>%
    mutate(ptime = (ptime - min(ptime)), 
           ptime = ptime / max(ptime))
  # invert if HSCs are > 0.5
  
  ggplot(dat.ptime.filt, aes(x = ptime)) + 
    geom_density() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  jcells <- subset(dat.glmpca.keep.filt, ctype.from.LL == "HSCs")$cell
  ptime.hscs <- median(subset(dat.ptime.filt, cell %in% jcells)$ptime)
  if (ptime.hscs > 0.5){
    # flip ptime
    print(paste("Flipping ptime for", jctype))
    dat.ptime.filt <- dat.ptime.filt %>%
      rowwise() %>%
      mutate(ptime = 1 - ptime)
  }
  
  # plot again
  dat.glmpca.keep.filt.annot <- dat.glmpca.keep.filt %>%
    left_join(., dat.ptime.filt)
  
  m <- ggplot(dat.glmpca.keep.filt.annot, aes_string(x = jdim1, y = jdim2, color = "ptime")) + 
    geom_point() + 
    scale_color_viridis_c(na.value = "grey85") + 
    theme_bw() + 
    ggtitle(jtitle) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  return(dat.glmpca.keep.filt.annot)
})

dev.off()
saveRDS(dat.traj.lst, file = outrds)









# Save new metadata ready for fitting -------------------------------------



