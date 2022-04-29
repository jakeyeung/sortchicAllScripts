# Jake Yeung
# Date of Creation: 2022-04-05
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/8-downstream_LL_celltyping_repressive_cleaned_update_ctypes.R
# 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


jmarks <- c("k4me1", "k9me3", "k27me3", "k4me3"); names(jmarks) <- jmarks
# jmarks <- c("k4me1", "k4me3"); names(jmarks) <- jmarks
jsuffix <- "dynamicbins"
# indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_repressive_cleaned_Tcells"
# indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_repressive_cleaned_update_ctypes"
indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA"
outdir <- indir

# jmark <- jmarks[[1]]


# Get colors --------------------------------------------------------------

dat.meta.colors <- fread("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_cleaned_up_update_k27me3_umap/metadata.k4me1.2022-03-01.txt") %>%
  dplyr::select(ctype.from.LL, colcode) 
dat.meta.colors <- dat.meta.colors[!duplicated(dat.meta.colors), ]

# replace monocyt ecolor with #d69941
colors.hash <- hash::hash(dat.meta.colors$ctype.from.LL, dat.meta.colors$colcode)
colors.hash[["Monocytes"]] < "d69941"


for (jmark in jmarks){
  print(jmark)
  
  outmeta <- file.path(outdir, paste0("metadata_reannotate_from_LLmat_fix_ctypes_by_batch2_", jsuffix, ".", jmark, ".", Sys.Date(), ".txt"))
  outpdf <- file.path(outdir, paste0("plots_reannotate_from_LLmat_fix_ctypes_by_batch2_", jsuffix, ".", jmark, ".", Sys.Date(), ".pdf"))
  pdf(outpdf, useDingbats = FALSE)
  
  inf <- file.path(indir, paste0("LLmat_by_batch_dynamicbins_", jmark, ".2022-04-06.rds"))
  
  LL.all <- readRDS(inf)
  
  LL.mat <- t(LL.all)
  
  # indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/repressive_cleaned"
  
  if (jmark %in% c("k4me1", "k4me3")){
    indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_cleaned_up_update_k27me3_umap"
    inf.meta <- file.path(indir.meta, paste0("metadata.", jmark, ".2022-03-01.txt"))
    dat.meta <- fread(inf.meta) %>%
      dplyr::select(-ctype.from.LL)
  } else if (jmark %in% c("k27me3", "k9me3")) {
    inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/repressive_cleaned_check_eryths/metadata_celltyping_", jmark, ".dynamicbins.2022-04-06.txt")
    dat.meta <- fread(inf.meta)
  }
  ctypes.indx <- colnames(LL.mat)
  
  LL.best <- apply(LL.mat, 1, function(jrow){
    ctypes.indx[which.max(jrow)]
  })
  
  LL.secondbest <- apply(LL.mat, 1, function(jrow){
    ctypes.indx[which.max(jrow[jrow!=max(jrow)])]
  })
  
  dat.LL.best <- data.frame(cell = names(LL.best), ctype.from.LL = LL.best, ctype.from.LL.secondbest = LL.secondbest, stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(ctype.from.LL = gsub("^Eryths.old$", "Eryths", ctype.from.LL), 
           ctype.from.LL = gsub("^Bcells.old$", "Bcells", ctype.from.LL),
           ctype.from.LL = gsub("^Granulocytes.old$", "Granulocytes", ctype.from.LL),
           ctype.from.LL.secondbest = gsub("^Eryths.old$", "Eryths", ctype.from.LL.secondbest), 
           ctype.from.LL.secondbest = gsub("^Bcells.old$", "Bcells", ctype.from.LL.secondbest),
           ctype.from.LL.secondbest = gsub("^Granulocytes.old$", "Granulocytes", ctype.from.LL.secondbest))
  
  dat.meta.reannotate <- dat.meta %>%
    left_join(., dat.LL.best)
  


# Adjusts -----------------------------------------------------------------
  
  # adjust basophils
  # dat.meta.reannotate <- dat.meta.reannotate %>%
  #   rowwise() %>%
  #   mutate(ctype.from.LL = ifelse(ctype.from.LL == "Basophils" & ctype != "Basophils", ctype.from.LL.secondbest, ctype.from.LL))
  
  # adjust progenitors
  dat.meta.reannotate <- dat.meta.reannotate %>%
    rowwise() %>%
    mutate(ctype.from.LL = ifelse(ctype == "MEP", ctype, ctype.from.LL),
           ctype.from.LL = ifelse(ctype == "CMP", ctype, ctype.from.LL), 
           ctype.from.LL = ifelse(ctype == "MPPs", ctype, ctype.from.LL), 
           ctype.from.LL = ifelse(ctype == "HSCs", ctype, ctype.from.LL))
  
  
  # adjust NKs and B cells
  dat.meta.reannotate <- dat.meta.reannotate %>%
    rowwise() %>%
    mutate(ctype.from.LL = ifelse(ctype.from.LL == "NKs" & ctype == "Bcells", ctype, ctype.from.LL))
  
  dat.meta.reannotate <- dat.meta.reannotate %>%
    rowwise() %>%
    mutate(ctype.from.LL = ifelse(ctype.from.LL == "Bcells" & ctype == "NKs", ctype, ctype.from.LL))
  
  # adjust basophils for H3K4me3
  if (jmark == "k4me3"){
    dat.meta.reannotate <- dat.meta.reannotate %>%
      rowwise() %>%
      mutate(ctype.from.LL.tmp = ifelse(ctype.from.LL == "Basophils" & ctype != "Basophils" & batch == "Old", ctype, NA),
             ctype.from.LL.tmp = gsub("^HSPCs$", "HSCs", ctype.from.LL.tmp), 
             ctype.from.LL.tmp = gsub("^DCs$", "Monocytes", ctype.from.LL.tmp), 
             ctype.from.LL = ifelse(!is.na(ctype.from.LL.tmp), ctype.from.LL.tmp, ctype.from.LL))
  } 
  
  # adjust pDCs  and NKs only in k27me3
  if (jmark == "k27me3"){
    dat.meta.reannotate <- dat.meta.reannotate %>%
      rowwise() %>%
      mutate(ctype.from.LL = ifelse(ctype.from.LL == "pDCs" & ctype != "pDCs", ctype.from.LL.secondbest, ctype.from.LL))
    
    # fix NK cells
    # dat.meta.reannotate <- dat.meta.reannotate %>%
    #   rowwise() %>%
    #   mutate(ctype.from.LL = ifelse(ctype.from.LL == "NKs" & umap2 < 30, ctype.from.LL.secondbest, ctype.from.LL))
    
    # dat.meta.reannotate <- dat.meta.reannotate %>%
    #   rowwise() %>%
    #   mutate(ctype.from.LL = ifelse(ctype.from.LL == "NKs" & ctype != "NKs", ctype.from.LL.secondbest, ctype.from.LL))
    
    # cells.check <- subset(dat.meta.reannotate, umap2 < 30 & ctype.from.LL == "NKs")$cell
    # table(subset(dat.meta.reannotate.test, cell %in% cells.check)$ctype.from.LL)
    # 
    # jcheck <- dat.meta.reannotate %>%
    #   rowwise() %>%
    #   mutate(highlight = umap2 < 30 & ctype.from.LL == "NKs", 
    #          ctype.check = ifelse(highlight, ctype.from.LL.secondbest, NA))
    # 
    # ggplot(jcheck, aes(x = umap1, y = umap2, color = ctype.check)) + 
    #   geom_point() + 
    #   facet_wrap(~ctype.check) + 
    #   theme_bw() + 
    #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
  }
  
  
  
  

# End adjusts -------------------------------------------------------------
  
  # dat.meta.reannotate <- subset(dat.meta.reannotate, select = -ctype.from.LL.secondbest)
  dat.meta.reannotate$colcode <- sapply(dat.meta.reannotate$ctype.from.LL, function(x) colors.hash[[x]])

  
  m1 <- ggplot(dat.meta.reannotate, aes(x = umap1, y = umap2, color = colcode)) +
    geom_point() + 
    theme_bw() + 
    ggtitle(paste(jmark, "New + Old")) + 
    facet_wrap(~ctype.from.LL) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m1)
  
  m2 <- ggplot(dat.meta.reannotate, aes(x = umap1, y = umap2, color = ctype)) +
    geom_point() + 
    theme_bw() + 
    ggtitle(paste(jmark, "New + Old")) + 
    facet_wrap(~ctype) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m2)
  
  m1.new <- ggplot(dat.meta.reannotate %>% filter(batch == "New"), aes(x = umap1, y = umap2, color = colcode)) +
    geom_point() + 
    theme_bw() + 
    ggtitle(paste(jmark, "New only")) + 
    facet_wrap(~ctype.from.LL) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m1.new)
  
  
  
  m2.new <- ggplot(dat.meta.reannotate %>% filter(batch == "New"), aes(x = umap1, y = umap2, color = ctype)) +
    geom_point() + 
    theme_bw() + 
    ggtitle(paste(jmark, "New only")) + 
    facet_wrap(~ctype) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m2.new)
  
  m1.old <- ggplot(dat.meta.reannotate %>% filter(batch == "Old"), aes(x = umap1, y = umap2, color = colcode)) +
    geom_point() + 
    theme_bw() + 
    ggtitle(paste(jmark, "Old only")) + 
    facet_wrap(~ctype.from.LL) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m1.old)
  
  m2.old <- ggplot(dat.meta.reannotate %>% filter(batch == "Old"), aes(x = umap1, y = umap2, color = ctype)) +
    geom_point() + 
    theme_bw() + 
    ggtitle(paste(jmark, "Old only")) + 
    facet_wrap(~ctype) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m2.old)
  
  
  
  
  # check DCs
  jctypes <- sort(unique(dat.meta.reannotate$ctype.from.LL))
  jbatches <- c("New", "Old")
  for (jctype in jctypes){
    print(jctype)
    for (jbatch in jbatches){
      print(jbatch)
      
      jsub <- dat.meta.reannotate %>% 
        mutate(highlight.ctype = ctype == jctype & batch == jbatch,
               highlight.ctype.LL = ctype.from.LL == jctype & batch == jbatch)
      m1.new2 <- ggplot(jsub %>%
                          arrange(highlight.ctype, highlight.ctype.LL), 
                        aes(x = umap1, y = umap2, color = highlight.ctype)) +
        geom_point() + 
        theme_bw() + 
        ggtitle(paste(jmark, jctype, jbatch)) + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
      m1.new2.LL <- ggplot(jsub %>%
                             arrange(highlight.ctype.LL), 
                           aes(x = umap1, y = umap2, color = highlight.ctype.LL)) +
        geom_point() + 
        theme_bw() + 
        ggtitle(paste(jmark, jctype, jbatch)) + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
      
      # ggplot(jsub %>% filter(highlight.ctype.LL), aes(x = umap1, y = umap2, color = ctype)) + 
      #   geom_point() + 
      #   facet_wrap(~ctype) + 
      #   theme_bw() + 
      #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      # 
      # ggplot(jsub %>% filter(highlight.ctype.LL), aes(x = umap1, y = umap2, color = ctype.from.LL.secondbest)) + 
      #   geom_point() + 
      #   facet_wrap(~ctype.from.LL.secondbest) + 
      #   theme_bw() + 
      #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      
      JFuncs::multiplot(m1.new2, m1.new2.LL, cols = 2)
      
    }
  }

  
   

  # Add colors --------------------------------------------------------------
  
  
  m1.final <- ggplot(dat.meta.reannotate, aes(x = umap1, y = umap2, color = colcode)) +
    geom_point() + 
    theme_bw() + 
    ggtitle(paste(jmark, "New + Old")) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m1.final)
  
  m1.final.batch <- ggplot(dat.meta.reannotate, aes(x = umap1, y = umap2, color = colcode)) +
    geom_point() + 
    theme_bw() + 
    facet_wrap(~batch) + 
    ggtitle(paste(jmark, "New + Old")) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m1.final.batch)
  
  
  # Write outputs -----------------------------------------------------------
  
  fwrite(dat.meta.reannotate, file = outmeta, sep = "\t")
  dev.off()
}



# 
# # annotate some of the celltypes
# 
# dat.meta.to.change <- subset(dat.meta.reannotate, ctype %in% c("AllCells", "IL7RLinNeg", "LinNeg", "LSK", "HSPCs"))  %>%
#   dplyr::mutate(ctype = ctype.from.LL)
# dat.meta.keep <- subset(dat.meta.reannotate, !ctype %in% c("AllCells", "IL7RLinNeg", "LinNeg", "LSK", "HSPCs"))
# 
# dat.meta.changed <- rbind(dat.meta.to.change, dat.meta.keep)
# 
# m2 <- ggplot(dat.meta.changed, aes(x = umap1, y = umap2, color = ctype)) +
#   geom_point() + 
#   theme_bw() + 
#   facet_wrap(~ctype) + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
# 
# JFuncs::multiplot(m1, m2, cols = 2)







