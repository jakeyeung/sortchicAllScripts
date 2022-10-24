# Jake Yeung
# Date of Creation: 2022-08-01
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/12-check_plate_effects.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("k4me1", "k4me3", "k27me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarksold) <- jmarks

# Load meta and plot  -----------------------------------------------------

indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/downstream_LDA_progs_only"

dat.meta.lst <- lapply(jmarks, function(jmark){
  # inf <- file.path(indir, paste0("LDA_downstream_celltypes_Giladi.", jmark, ".pdf"))
  inf <- file.path(indir, paste0("dat_umap_progs_only.", jmark, ".2022-07-29.txt"))
  fread(inf) %>%
    rowwise() %>%
    mutate(plate = scchicFuncs::ClipLast(cell, jsep = "_"), 
           experi = scchicFuncs::ClipLast(plate, jsep = "-"))
})

outpdf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/downstream_LDA_progs_only/check_plate_effects.", Sys.Date(), ".pdf")
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
pdf(outpdf, useDingbats = FALSE)
m.batch <- lapply(jmarks, function(jmark){
  jdat <- dat.meta.lst[[jmark]]
  ggplot(jdat, aes(x = umap1, y = umap2, color = experi)) + 
    geom_point() + 
    ggtitle(jmark) +
    scale_color_manual(values = cbPalette) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
print(m.batch)
dev.off()

