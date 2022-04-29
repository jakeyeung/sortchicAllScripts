# Jake Yeung
# Date of Creation: 2022-04-20
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/24-check_PCA_of_LDA_downstream_other_marks.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)

jmarksnew <- c("k4me1", "k4me3", "k9me3"); names(jmarksnew) <- jmarksnew

# Load LDAs ---------------------------------------------------------------

indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths2"
outs <- lapply(jmarksnew, function(jmark){
  
  inf.meta <- file.path(indir.meta, paste0("metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.", jmark, ".txt"))
  if (jmark == "k4me1"){
    inf.ldaout <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.k4me1.2022-01-28/ldaOut.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.k4me1.2022-01-28.Robj"
  } else if (jmark == "k4me3"){
    inf.ldaout <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k4me3_cleaned/lda_outputs.count_mat_cleaned_no3no7_dynbins_allcells.k4me3.2022-04-12/ldaOut.count_mat_cleaned_no3no7_dynbins_allcells.k4me3.2022-04-12.Robj"
  } else if (jmark == "k27me3"){
    inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k27me3_clean_eryths/lda_outputs.count_mat_merged_with_old_dynbins.k27me3.2022-04-15/ldaOut.count_mat_merged_with_old_dynbins.", jmark, ".2022-04-15.Robj")
  } else if (jmark == "k9me3"){
    inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_repressive_cleaned_from_jupyter/ldaAnalysis_fripfilt_varfilt_binfilt/lda_outputs.count_mat_cleaned_dynbins.", jmark, ".2022-02-16/ldaOut.count_mat_cleaned_dynbins.", jmark, ".2022-02-16.Robj")
  }
  
  load(inf.ldaout, v=T)
  dat.meta <- fread(inf.meta)
  
  tm <- posterior(out.lda)
  
  return(list(tm = tm, dat.meta = dat.meta))
})



# Run PCAs ----------------------------------------------------------------

dat.impute.log.lst <- lapply(outs, function(out){
  t(log2(out$tm$topics %*% out$tm$terms))
})

dat.meta.lst <- lapply(outs, function(out){
  out$dat.meta
})


# Save PCAs ---------------------------------------------------------------

library(irlba)
ntopics <- 30
dat.pca.lst <- parallel::mclapply(dat.impute.log.lst, function(logodds){
  # remove mean and SVD
  logodds.centered <- t(scale(t(logodds), center = TRUE, scale = TRUE))
  # logodds.centered.check <- sweep(logodds, MARGIN = 1, STATS = rowMeans(logodds), FUN = "-")
  # logodds.pca <- prcomp(t(logodds.centered), center = FALSE, scale. = FALSE, rank. = ntopics)
  # logodds.pca <- prcomp(t(logodds.centered), center = FALSE, scale. = FALSE, rank. = ntopics)
  logodds.pca <- irlba(A = t(logodds.centered), nv = ntopics, scale = FALSE, center = FALSE)
  
  rownames(logodds.pca$u) <- colnames(logodds)
  rownames(logodds.pca$v) <- rownames(logodds)
  
  return(logodds.pca)
}, mc.cores = length(jmarksnew))


outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections"
for (jmark in jmarksnew){
  outftmp <- file.path(outdir, paste0("pca_output.irlba.scaled.", jmark, ".2022-04-20.rds"))
  dat.pca <- dat.pca.lst[[jmark]]
  saveRDS(dat.pca, file = outftmp)
}

# Visualize PCAs ----------------------------------------------------------


PoV.lst <- lapply(dat.pca.lst, function(logodds.pca){
  signif(logodds.pca$d^2/sum(logodds.pca$d^2), digits = 2)
})


dat.pca.annot.lst <- lapply(jmarksnew, function(jmark){
  # U.init <- dat.pca.lst[[jmark]]$x
  U.init <- dat.pca.lst[[jmark]]$u
  rownames(U.init) <- colnames(dat.impute.log.lst[[jmark]])
  dat.pca <- data.frame(cell = rownames(U.init), pc1 = U.init[, 1], pc2 =U.init[, 2], pc3 = U.init[, 3], pc4 = U.init[, 4], stringsAsFactors = FALSE) %>%
    left_join(., dat.meta.lst[[jmark]])
})

dat.meta.colors <- subset(dat.meta.lst[[1]], select = c(colcode, ctype.from.LL))
dat.meta.colors <- dat.meta.colors[!duplicated(dat.meta.colors), ]

m.lst <- lapply(jmarksnew, function(jmark){
  jdat <- dat.pca.annot.lst[[jmark]] %>%
    mutate(ctype.order = ifelse(ctype.from.LL == "HSCs", "zzzHSCs", ctype.from.LL)) %>%
    mutate(ctype.order = ifelse(ctype.order == "LT", "zzLT", ctype.order)) %>%
    mutate(ctype.order = ifelse(ctype.order == "ST", "zST", ctype.order)) %>%
    arrange(ctype.order)
  m <- ggplot(jdat, aes(x = pc1, y = pc2, color = colcode)) + 
    geom_point()  + 
    theme_bw() + 
    ggtitle(jmark) + 
    xlab(paste0("PC1 (", PoV.lst[[jmark]][[1]], ")")) + 
    ylab(paste0("PC2 (", PoV.lst[[jmark]][[2]], ")")) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m)
})

outpdf <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections"
pdf(file.path(outpdf, paste0("pca_other_marks.", Sys.Date(), ".pdf")), useDingbats = FALSE)
print(m.lst)
dev.off()

