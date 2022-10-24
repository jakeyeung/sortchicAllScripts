# Jake Yeung
# Date of Creation: 2022-07-27
# File: ~/projects/scchic/scripts/revision_scripts/revisions2/3-Liger_downstreams.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(rliger)

jmarks <- c("k4me1", "k4me3", "k27me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarksold) <- jmarks


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

dat.meta.merge <- lapply(jmarks, function(jmark){
  jdat <- dat.meta.lst[[jmark]]
  subset(jdat, select = c(cell, ctype.from.LL, colcode)) %>%
    mutate(mark = jmark)
}) %>%
  bind_rows()


# Get marker genes --------------------------------------------------------

dat.markers <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/genesets/dat_genes_sub_join_from_orig.rds")
print(unique(dat.markers$jset))

genes.keep <- unique(dat.markers$gene)



# Ligers ------------------------------------------------------------------




indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/liger_outputs"

jmarkref <- "k4me3"
# jmarkref <- "k4me1"
jmarkothers <- jmarks[jmarks != jmarkref]

outpdf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/integration_data/liger_outputs/plots/", jmarkref, ".", Sys.Date(), ".pdf")
pdf(outpdf, useDingbats = FALSE)

ligerobj.lst <- lapply(jmarkothers, function(jmark){
  print(jmark)
  ligerobj <- readRDS(file.path(indir, paste0("liger_", jmarkref, "_and_", jmark, "_UMAP.2022-07-27.rds")))
  dat.umap.liger <- data.frame(cell = rownames(ligerobj@tsne.coords), ligerobj@tsne.coords, stringsAsFactors = FALSE) %>%
    left_join(., dat.meta.merge)
})

m.lst <- lapply(jmarkothers, function(jmark){
  dat.meta.liger <- ligerobj.lst[[jmark]]
  if (jmarkref == "k4me1"){
    dat.meta.liger %>% filter(X1 > -14)
  }
  ggplot(dat.meta.liger, aes(x = X1, y = X2, color = colcode)) + 
    geom_point() + 
    ggtitle(paste(jmarkref, "and", jmark)) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode, guide = "legend") +
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
print(m.lst)

m.facet.lst <- lapply(jmarkothers, function(jmark){
  dat.meta.liger <- ligerobj.lst[[jmark]]
  if (jmarkref == "k4me1"){
    dat.meta.liger %>% filter(X1 > -14)
  }
  ggplot(dat.meta.liger, aes(x = X1, y = X2, color = colcode)) + 
    geom_point() + 
    facet_wrap(~ctype.from.LL) + 
    ggtitle(paste(jmarkref, "and", jmark)) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode, guide = "legend") +
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
print(m.facet.lst)


m.mark.lst <- lapply(jmarkothers, function(jmark){
  dat.meta.liger <- ligerobj.lst[[jmark]]
  if (jmarkref == "k4me1"){
    dat.meta.liger %>% filter(X1 > -14)
  }
  ggplot(dat.meta.liger, aes(x = X1, y = X2, color = mark)) + 
    geom_point(alpha = 0.5) + 
    ggtitle(paste(jmarkref, "and", jmark)) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
print(m.mark.lst)


m.mark.lst <- lapply(jmarkothers, function(jmark){
  dat.meta.liger <- ligerobj.lst[[jmark]]
  if (jmarkref == "k4me1"){
    dat.meta.liger %>% filter(X1 > -14)
  }
  ggplot(dat.meta.liger, aes(x = X1, y = X2, color = mark)) + 
    geom_point(alpha = 0.5) + 
    ggtitle(paste(jmarkref, "and", jmark)) + 
    facet_wrap(~mark) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
print(m.mark.lst)

dev.off()