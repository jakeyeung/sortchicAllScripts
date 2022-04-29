# Jake Yeung
# Date of Creation: 2022-04-14
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/9-pseudotime_for_different_celltypes_k4me1_cleaned.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmark <- "k4me1"


# Load metadata -----------------------------------------------------------

# inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping/metadata_reannotate_from_LLmat_dynamicbins.", jmark, ".2022-02-02.txt")
# inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping/metadata_reannotate_from_LLmat_dynamicbins.", jmark, ".2022-02-03.txt")
inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths/metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.", jmark, ".2022-04-13.txt")
dat.meta.LL <- fread(inf.meta)


# Load raw counts ---------------------------------------------------------

inf.glmpca <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/glmpca_outputs/varfilt/glmpca.", jmark, ".bincutoff_0.binskeep_0.platename_batch.szname_none.niter_1000.RData")
load(inf.glmpca, v=T)
count.mat <- glm.inits$Y.filt

dat.glmpca <- data.frame(cell = rownames(glm.out$factors), glm.out$factors, stringsAsFactors = FALSE) %>%
  left_join(., dat.meta.LL)


# Create different trajectories -------------------------------------------

hpcs <- c("HSCs", "LT", "ST", "MPPs")
# granu tracj
ctypes.granus <- c(hpcs, "CMP", "GMP", "Granulocytes")
dims.granus <- c("dim1", "dim2")

ctypes.bcells <- c(hpcs, "Bcells")
dims.bcells <- c("dim7", "dim8")

ctypes.eryths <- c(hpcs, "MEP", "Eryths")
dims.eryths <- c("dim2", "dim3")

ctypes.pdcs <- c(hpcs, "pDCs")
dims.pdcs <- c("dim8", "dim9")

ctypes.dcs <- c(hpcs, "CMP", "GMP", "DCs")
dims.dcs <- c("dim1", "dim2")

ctypes.basos <- c(hpcs, "CMP", "GMP", "Basophils")
dims.basos <- c("dim1", "dim2")

ctypes.nks <- c(hpcs, "NKs")
dims.nks <- c("dim4", "dim5")

ctypes.monos <- c(hpcs, "CMP", "GMP", "Monocytes")
dims.monos <- c("dim1", "dim2")

ctypes.lst <- list(ctypes.granus, ctypes.bcells, ctypes.eryths, ctypes.pdcs, ctypes.dcs, ctypes.basos, ctypes.nks, ctypes.monos)
names(ctypes.lst) <- sapply(ctypes.lst, function(x) x[length(x)])

ctypes.names <- names(ctypes.lst); names(ctypes.names) <- ctypes.names

dims.lst  <- list(dims.granus, dims.bcells, dims.eryths, dims.pdcs, dims.dcs, dims.basos, dims.nks, dims.monos)
names(dims.lst) <- ctypes.names

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/cleaned"
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
  
  m <- ggplot(dat.glmpca.keep, aes_string(x = jdim1, y = jdim2, color = "is.ctype")) + 
    geom_point() + 
    theme_bw() + 
    # facet_wrap(~ctype.from.LL) + 
    ggtitle(jtitle) + 
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
    ungroup() %>%
    mutate(ctype.ordered = factor(ctype.from.LL, levels = ctypes.vec))
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



