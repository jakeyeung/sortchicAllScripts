# Jake Yeung
# Date of Creation: 2022-01-28
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/3-check_celltype_labeling.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

BatchColumn2Ctype <- function(Batch, Column, Row = NA){
  
  if (Batch == "SL1"){
    if (Column >= 1 & Column <= 4){
      ctype <- "Tcells"
    } else if (Column >= 5 & Column <= 8){
      ctype <- "Bcells"
    } else if (Column >= 9 & Column <= 12){
      ctype <- "NKs"
    } else if (Column >= 13 & Column <= 16){
      ctype <- "Eryths"
    } else if (Column >= 17 & Column <= 24)
      ctype <- "AllCells"
  }
  
  if (Batch == "SL2"){
    if (Column >= 1 & Column <= 4){
      ctype <- "Granulocytes"
    } else if (Column >= 5 & Column <= 8){
      ctype <- "Monocytes"
    } else if (Column >= 9 & Column <= 12){
      ctype <- "DCs"
    } else if (Column >= 13 & Column <= 24){
      ctype <- "AllCells"
    }
  }
  
  if (Batch == "SL3"){
    if (Column == 1){
      ctype <- "HSCs"
    } else if (Column == 2){
      ctype <- "MPPs"
    } else if (Column == 3){
      ctype <- "LT"
    } else if (Column == 4){
      ctype <- "ST"
    } else if (Column >= 5 & Column <= 24){
      ctype <- "LSK"
    } else {
      warning("Error: unknown column:", Column)
    }
  }
  
  if (Batch == "SL4"){
    if (Column == 1){
      ctype <- "GMP"
    } else if (Column == 2){
      ctype <- "CMP"
    } else if (Column == 3){
      ctype <- "MEP"
    } else if (Column >= 4 & Column <= 24){
      ctype <- "LSK"
    }
  }
  
  if (Batch == "SL5"){
    if (Column == 1 & Row >= 1 & Row <= 8){
      ctype <- "pDCs"
    } else if (Column == 1 & Row >= 9 & Row <= 16){
      ctype <- "IL7RLinNeg"
    } else if (Column >= 2 & Column <= 12){
      ctype <- "LinNeg"
    } else if (Column >= 13 & Column <= 24){
      ctype <- "LSK"
    }
  }
  
  return(ctype)
}

GetPlateCoord2 <- function(cell, platecols = 24, is.zero.base = FALSE){
  # cell0 -> 1,1
  indx <- as.numeric(strsplit(cell, "cell")[[1]][[2]]) 
  if (is.zero.base){
    indx <- indx + 1
  }
  jcol <- indx %% platecols
  jcol <- ifelse(jcol == 0, platecols, jcol)
  jrow <- ceiling(indx / platecols)
  return(c(jrow, jcol))
}


indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/metadata"

jsuffix1 <- "_dynamic_bins"
jsuffix2 <- "_merged_with_old"

jmarksnew <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarksnew) <- jmarksnew
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksold) <- jmarksnew

# jmarkold <- "H3K4me1"
# jmark <- "k4me1"
# print(jmark)

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping"
for (jmark in jmarksnew){
  print(jmark)
  
  outpdf <- file.path(outdir, paste0("celltyping_", jmark, jsuffix1, jsuffix2, ".", Sys.Date(), ".pdf"))
  outmeta <- file.path(outdir, paste0("metadata_celltyping_", jmark, jsuffix1, jsuffix2, ".", Sys.Date(), ".txt"))
  
  pdf(outpdf, useDingbats = FALSE)
  
  
  jmarkold <- jmarksold[[jmark]]
  
  
  indir <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_bugfixed", jsuffix1)
  inf <- file.path(indir, paste0("ldaAnalysis_fripfilt_BM_", jmark, "/lda_outputs.count_mat", jsuffix2, ".", jmark, ".2022-01-26/ldaOut.count_mat", jsuffix2, ".", jmark, ".2022-01-26.Robj"))
  assertthat::assert_that(file.exists(inf))
  
  
  inf.meta.old <- file.path(indir.meta, paste0("metadata_batch_corrected.arranged_by_lineage.", jmarkold, ".2021-01-02.txt"))
  inf.meta.new <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata/metadata_plate_experi_batch.", jmark, ".2022-01-27.txt")
  
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
    ggtitle(paste(jmark, jsuffix1, jsuffix2)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.meta.new, aes(x = umap1, y = umap2, color = batch)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(paste(jmark, jsuffix1, jsuffix2)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
  
  m <- ggplot(dat.meta.new %>% filter(batch == "Old"), aes(x = umap1, y = umap2, color = ctype)) + 
    geom_point() + 
    facet_wrap(~ctype) + 
    theme_bw() + 
    ggtitle(paste(jmark, jsuffix1, jsuffix2)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
  m <- ggplot(dat.meta.new %>% filter(batch == "New"), aes(x = umap1, y = umap2, color = ctype)) + 
    geom_point() + 
    facet_wrap(~ctype) + 
    theme_bw() + 
    ggtitle(paste(jmark, jsuffix1, jsuffix2)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  ctypes.highlight <- sort(unique(subset(dat.meta.new, batch == "New")$ctype))
  
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



