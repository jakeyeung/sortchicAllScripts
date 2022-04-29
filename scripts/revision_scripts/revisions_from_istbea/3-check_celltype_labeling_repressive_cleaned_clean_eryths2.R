# Jake Yeung
# Date of Creation: 2022-04-13
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/3-check_celltype_labeling_repressive_cleaned_keep_eryths.R
# 



rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

indir.meta.old <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/metadata"

jsuffix <- "dynamicbins"

jmarksnew <- c("k27me3"); names(jmarksnew) <- jmarksnew
jmarksold <- c("H3K27me3"); names(jmarksold) <- jmarksnew

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/repressive_cleaned_clean_eryths2"
dir.create(outdir)

jdate <- Sys.Date() 

for (jmark in jmarksnew){
  print(jmark)
  
  outpdf <- file.path(outdir, paste0("celltyping_", jmark, ".", jsuffix, ".", jdate, ".pdf"))
  outmeta <- file.path(outdir, paste0("metadata_celltyping_", jmark, ".", jsuffix, ".", jdate, ".txt"))
  
  pdf(outpdf, useDingbats = FALSE)
  
  
  jmarkold <- jmarksold[[jmark]]
  
  
  inf.meta.old <- file.path(indir.meta.old, paste0("metadata_batch_corrected.arranged_by_lineage.", jmarkold, ".2021-01-02.txt"))
  # inf.meta.new <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata/repressive_cleaned_keep_eryths/metadata_plate_experi_batch.merged_with_old_dynbins.dynamicbins.", jmark, ".2022-04-13.txt")
  inf.meta.new <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata/repressive_cleaned_clean_eryths2/metadata_plate_experi_batch.merged_with_old_dynbins.dynamicbins.", jmark, ".2022-04-16.txt")
  assertthat::assert_that(file.exists(inf.meta.new))
  
  dat.meta.old <- fread(inf.meta.old)
  dat.meta.new <- fread(inf.meta.new) %>%
    rowwise() %>%
    mutate(indx = as.numeric(strsplit(cell, split = "_")[[1]][[2]]),
           cellindx = paste0("cell", indx), 
           platerow = GetPlateCoord2(cell = cellindx, platecols = 24, is.zero.base = FALSE)[[1]], 
           platecol = GetPlateCoord2(cell = cellindx, platecols = 24, is.zero.base = FALSE)[[2]], 
           Batch = ifelse(batch == "New", strsplit(experi, split = "-")[[1]][[4]], "Old"))
  
  cell2ctype <- hash::hash(dat.meta.old$cell, dat.meta.old$cluster)
  
  dat.meta.new <- dat.meta.new %>%
    rowwise() %>%
    mutate(ctype = ifelse(Batch == "Old", scchicFuncs::AssignHash(x = cell, jhash = cell2ctype, null.fill = "Error"), BatchColumn2Ctype(Batch = Batch, Column = platecol, Row = platerow)))
  
  m <- ggplot(dat.meta.new, aes(x = umap1, y = umap2, color = Batch)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(paste(jmark, jsuffix)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.meta.new, aes(x = umap1, y = umap2, color = batch)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(paste(jmark, jsuffix)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
  
  m <- ggplot(dat.meta.new %>% filter(batch == "Old"), aes(x = umap1, y = umap2, color = ctype)) + 
    geom_point() + 
    facet_wrap(~ctype) + 
    theme_bw() + 
    ggtitle(paste(jmark, jsuffix)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
  m <- ggplot(dat.meta.new %>% filter(batch == "New"), aes(x = umap1, y = umap2, color = ctype)) + 
    geom_point() + 
    facet_wrap(~ctype) + 
    theme_bw() + 
    ggtitle(paste(jmark, jsuffix)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  # ctypes.highlight <- sort(unique(subset(dat.meta.new, batch == "New")$ctype))
  ctypes.highlight <- sort(unique(dat.meta.new$ctype))
  
  print(ctypes.highlight)
  
  for (ctype.highlight in ctypes.highlight){
    print(ctype.highlight)
    
    dat.meta.new.tmp <- dat.meta.new %>% 
      rowwise() %>%
      mutate(ctype = ctype == ctype.highlight)
    
    subset(dat.meta.new.tmp, batch == "New")
    
    m <- ggplot(dat.meta.new.tmp %>% arrange(ctype), 
           aes(x = umap1, y = umap2, color = ctype)) + 
      geom_point() + 
      facet_wrap(~batch) + 
      theme_bw() + 
      ggtitle(paste(ctype.highlight)) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
    
    m <- ggplot(dat.meta.new %>%
             filter(ctype == ctype.highlight) %>%
             group_by(louvain) %>%
             summarise(ncells = length(cell)) %>%
             ungroup() %>%
             mutate(nfrac = ncells / sum(ncells)), 
           aes(x = as.character(louvain), y = nfrac)) + 
      geom_col() + 
      ggtitle(paste(ctype.highlight)) + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
    
    m <- ggplot(dat.meta.new.tmp %>% arrange(ctype), 
           aes(x = umap1, y = umap2, color = ctype)) + 
      geom_point() + 
      facet_wrap(~louvain) + 
      ggtitle(paste(ctype.highlight)) + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
  }
  fwrite(dat.meta.new, file = outmeta, sep = "\t")
  dev.off()
}



