# Jake Yeung
# Date of Creation: 2022-04-05
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/4-glmpca_downstream_repressive_cleaned.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(glmpca)


library(hash)
library(igraph)
library(umap)
library(scchicFuncs)


# Assign ctype ------------------------------------------------------------


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



# Settings ----------------------------------------------------------------



jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 30
jsettings[["min_dist"]] <- 0.1
jsettings[["random_state"]] <- 123
jsettings[["spread"]] <- 7

# jmark <- "k27me3"
outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/glmpca_outputs/repressive_cleaned"
dir.create(outdir)
jmarks <- c("k27me3", "k9me3"); names(jmarks) <- jmarks

for (jmark in jmarks){
  
  
  # Load glmpca  ------------------------------------------------------------
  
  # inf.glmpca <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/glmpca_outputs/glmpca.", jmark, ".bincutoff_0.binskeep_0.platename_batch.szname_none.niter_1000.RData")
  inf.glmpca <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/glmpca_outputs/repressive_cleaned/glmpca.", jmark, ".bincutoff_0.binskeep_0.platename_batch.szname_none.niter_1000.RData")
  load(inf.glmpca, v=T)
  
  bname <- basename(inf.glmpca)
  pdfname <- paste0(ClipLast(x = bname, jsep = "\\.", jsep.out = "."), ".pdf")
  outpdf <- file.path(outdir, pdfname)
  
  
  # Load meta  --------------------------------------------------------------
  
  # inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata/metadata_plate_experi_batch._dynamic_bins_merged_with_old.", jmark, ".2022-01-28.txt")
  # inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/metadata_celltyping_", jmark, "_dynamic_bins_merged_with_old.2022-01-28.txt")
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes/metadata_reannotate_from_LLmat_fix_pDCs_k27me3only_dynamicbins.", jmark, ".2022-04-05.txt")
  dat.meta <- fread(inf.meta)
  
  
  # Show UMAP  --------------------------------------------------------------
  
  dat.umap <- DoUmapAndLouvain(glm.out$factors, jsettings = jsettings)
  
  dat.meta.filt <- subset(dat.meta, select = c(cell, ctype, batch, Batch, plate, colcode, ctype.from.LL))
  dat.umap.annot <- left_join(dat.umap, dat.meta.filt)
    # rowwise() %>%
    # mutate(indx = as.numeric(strsplit(cell, split = "_")[[1]][[2]]),
    #        cellindx = paste0("cell", indx),
    #        experi = ClipLast(cell, jsep = "_"), 
    #        platerow = GetPlateCoord2(cell = cellindx, platecols = 24, is.zero.base = FALSE)[[1]],
    #        platecol = GetPlateCoord2(cell = cellindx, platecols = 24, is.zero.base = FALSE)[[2]],
    #        batch = ifelse(grepl("PZ-sortChIC-BM-SL", experi), "New", "Old"), 
    #        Batch = ifelse(batch == "New", strsplit(experi, split = "-")[[1]][[4]], "Old"),
    #        ctype = ifelse(Batch == "Old", ctype, BatchColumn2Ctype(Batch = Batch, Column = platecol, Row = platerow)))
  
  
  dat.meta.colors <- subset(dat.umap.annot, select = c(ctype.from.LL, colcode))
  dat.meta.colors <- dat.meta.colors[!duplicated(dat.meta.colors), ]
  
  m.umap.batch1 <- ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = batch)) + 
    geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
  m.umap.colcode <- ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = colcode)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
  m.umap.batch2 <- ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = batch)) + 
    geom_point() + 
    facet_wrap(~ctype) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
  
  # Show individual glmpcas -------------------------------------------------
  
  glm.factors <- data.frame(cell = rownames(glm.out$factors), glm.out$factors, stringsAsFactors = FALSE) %>%
    left_join(., dat.meta.filt)
  
  dims <- seq(ncol(glm.out$factors) - 1)
  # jfac1 <- "dim3"
  # jfac2 <- "dim4"
  
  pdf(outpdf, useDingbats = FALSE)
  print(m.umap.batch1)
  print(m.umap.colcode)
  print(m.umap.batch2)
  for (d1 in dims){
    d2 <- d1 + 1
    jfac1 <- paste0("dim", d1)
    jfac2 <- paste0("dim", d2)
    
    m1 <- ggplot(glm.factors, aes_string(x = jfac1, y = jfac2, color = "colcode")) + 
      geom_point() + 
      theme_bw() + 
      scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                            guide = "legend") + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    
    m2 <- ggplot(glm.factors, aes_string(x = jfac1, y = jfac2, color = "batch")) + 
      geom_point() + 
      theme_bw() + 
      facet_wrap(~ctype.from.LL) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
    
    print(m1)
    print(m2)
  }
  dev.off()
  
  
  
}
