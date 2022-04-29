# Jake Yeung
# Date of Creation: 2022-04-20
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/25-last_main_fig.R
# Assemble for last main fig

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks

# Load metas --------------------------------------------------------------

indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths2"

dat.meta.lst <- lapply(jmarks, function(jmark){
  inf.meta <- file.path(indir.meta, paste0("metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.", jmark, ".txt"))
  fread(inf.meta)
})


# Load UMAPs --------------------------------------------------------------

# for K27me3 we have a new UMAP 
inf.meta.k27me3.bc <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/dat_umap_annot_batch_corrected.k27me3.2022-04-21.txt"
dat.meta.k27me3.bc <- fread(inf.meta.k27me3.bc)
dat.meta.lst$k27me3 <- dat.meta.k27me3.bc

dat.meta.colors <- subset(dat.meta.k27me3.bc, select = c(ctype.from.LL, colcode))
dat.meta.colors <- dat.meta.colors[!duplicated(dat.meta.colors), ]

# scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode, guide = "legend")

m.lst <- lapply(jmarks, function(jmark){
  ggplot(dat.meta.lst[[jmark]], aes(x = umap1, y = umap2, color = colcode)) + 
    geom_point() + 
    theme_minimal() + 
    ggtitle(paste(jmark, "Ncells:", nrow(dat.meta.lst[[jmark]]))) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})


m.batch.lst <- lapply(jmarks, function(jmark){
  ggplot(dat.meta.lst[[jmark]], aes(x = umap1, y = umap2, color = colcode)) + 
    geom_point() + 
    theme_minimal() + 
    facet_wrap(~batch) + 
    ggtitle(paste(jmark, "Ncells:", nrow(dat.meta.lst[[jmark]]))) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})


# Load PCAs ---------------------------------------------------------------

indir.pca <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections"

pca.out.lst <- lapply(jmarks, function(jmark){
  pca.out <- readRDS(file.path(indir.pca, paste0("pca_output.", jmark, ".rds")))
})

PoV.lst <- lapply(jmarks, function(jmark){
  PoV <- signif(pca.out.lst[[jmark]]$d^2/sum(pca.out.lst[[jmark]]$d^2), digits = 2)
})

dat.pca.lst <- lapply(jmarks, function(jmark){
  dat.tmp <- pca.out.lst[[jmark]]$u
  colnames(dat.tmp) <- paste("PC", seq(ncol(dat.tmp)), sep = "")
  dat.pca <- data.frame(cell = rownames(dat.tmp), dat.tmp, stringsAsFactors = FALSE) %>%
    left_join(., dat.meta.lst[[jmark]], by = "cell")
  if (jmark %in% c("k4me1", "k4me3")){
    dat.pca$PC1 <- dat.pca$PC1 * -1
  }
  return(dat.pca)
})

m.pca.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.pca.lst[[jmark]] %>%
    mutate(ctype.order = ifelse(ctype.from.LL == "HSCs", "zzzHSCs", ctype.from.LL)) %>%
    mutate(ctype.order = ifelse(ctype.order == "LT", "zzLT", ctype.order)) %>%
    mutate(ctype.order = ifelse(ctype.order == "ST", "zST", ctype.order)) %>%
    arrange(ctype.order)
  ggplot(jdat, aes(x = PC1, y = PC2, color = colcode)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(paste(jmark, "Ncells:", nrow(dat.meta.lst[[jmark]]))) + 
    xlab(paste0("PC1 (", PoV.lst[[jmark]][[1]], ")")) + 
    ylab(paste0("PC2 (", PoV.lst[[jmark]][[2]], ")")) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    
})


m.pca.batch.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.pca.lst[[jmark]] %>%
    mutate(ctype.order = ifelse(ctype.from.LL == "HSCs", "zzzHSCs", ctype.from.LL)) %>%
    mutate(ctype.order = ifelse(ctype.order == "LT", "zzLT", ctype.order)) %>%
    mutate(ctype.order = ifelse(ctype.order == "ST", "zST", ctype.order)) %>%
    arrange(ctype.order)
  ggplot(jdat, aes(x = PC1, y = PC2, color = colcode)) + 
    geom_point() + 
    theme_bw() + 
    facet_wrap(~batch) + 
    ggtitle(paste(jmark, "Ncells:", nrow(dat.meta.lst[[jmark]]))) + 
    xlab(paste0("PC1 (", PoV.lst[[jmark]][[1]], ")")) + 
    ylab(paste0("PC2 (", PoV.lst[[jmark]][[2]], ")")) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
})

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umaps_pcas_with_batch_corrections"
pdf(file.path(outdir, paste0("umaps_pcas_all_marks_k27me3_batch_corrected.", Sys.Date(), ".pdf")))

print(m.lst)
print(m.batch.lst)

print(m.pca.lst)
print(m.pca.batch.lst)

dev.off()


# Save objects ------------------------------------------------------------

for (jmark in jmarks){
  outmeta.tmp <- file.path(outdir, paste0("umap_metadata_primetime.", jmark, ".", Sys.Date(), ".txt"))
  outpca.tmp <- file.path(outdir, paste0("pca_metadata_primetime.", jmark, ".", Sys.Date(), ".txt"))
  
  fwrite(dat.meta.lst[[jmark]], file = outmeta.tmp)
  fwrite(dat.pca.lst[[jmark]], file = outpca.tmp)
  
}

