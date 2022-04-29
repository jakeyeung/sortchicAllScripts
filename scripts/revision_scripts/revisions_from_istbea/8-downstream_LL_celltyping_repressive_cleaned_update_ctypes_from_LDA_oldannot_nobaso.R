# Jake Yeung
# Date of Creation: 2022-04-10
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/8-downstream_LL_celltyping_repressive_cleaned_update_ctypes_from_LDA_oldannot.R
# 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Clean stuff up ----------------------------------------------------------

RemoveOldFromCtype <- function(ctype.from.LL){
  ctype.from.LL <- gsub("^Eryths.old$", "Eryths", ctype.from.LL)
  ctype.from.LL <- gsub("^Bcells.old$", "Bcells", ctype.from.LL)
  ctype.from.LL <- gsub("^Granulocytes.old$", "Granulocytes", ctype.from.LL)
  return(ctype.from.LL)
}


# marks -------------------------------------------------------------------

jmarks <- c("k27me3", "k9me3", "k4me1", "k4me3"); names(jmarks) <- jmarks

jmarks.active <- c("k4me1", "k4me3"); names(jmarks.active)
jmarks.repress <- c("k27me3", "k9me3"); names(jmarks.repress)

# jmarks <- c("k4me1", "k4me3"); names(jmarks) <- jmarks
jsuffix <- "dynamicbins"
indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_oldannot_nobasos"
outdir <- indir

# jmark <- jmarks[[1]]

prog.ctypes <- c("CMP", "MEP", "HSCs", "MPPs"); names(prog.ctypes) <- prog.ctypes

# Get colors --------------------------------------------------------------

dat.meta.colors <- fread("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_cleaned_up_update_k27me3_umap/metadata.k4me1.2022-03-01.txt") %>%
  dplyr::select(ctype.from.LL, colcode) 
dat.meta.colors <- dat.meta.colors[!duplicated(dat.meta.colors), ]

# replace monocyt ecolor with #d69941
colors.hash <- hash::hash(dat.meta.colors$ctype.from.LL, dat.meta.colors$colcode)
colors.hash[["Monocytes"]] < "d69941"


for (jmark in jmarks){
  print(jmark)
  
  if (jmark %in% jmarks.active){
    indir.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_", jsuffix, ".out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28")
    inf.ldaout <- file.path(indir.ldaout, paste0("ldaOut.count_mat_var_filt_", jsuffix, ".out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28.Robj"))
    load(inf.ldaout, v=T)
  } else if (jmark %in% jmarks.repress){
    inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_repressive_cleaned_from_jupyter/ldaAnalysis_fripfilt_varfilt_binfilt/lda_outputs.count_mat_cleaned_dynbins.", jmark, ".2022-02-16/ldaOut.count_mat_cleaned_dynbins.", jmark, ".2022-02-16.Robj")
    load(inf.ldaout, v=T)
  }
  
  
  outmeta <- file.path(outdir, paste0("metadata_reannotate_from_LLmat_fix_ctypes_by_batch_nobaso_", jsuffix, ".", jmark, ".", Sys.Date(), ".txt"))
  outpdf <- file.path(outdir, paste0("plots_reannotate_from_LLmat_fix_ctypes_by_batch_nobaso_", jsuffix, ".", jmark, ".", Sys.Date(), ".pdf"))
  pdf(outpdf, useDingbats = FALSE)
  
  inf <- file.path(indir, paste0("LLmat_by_batch_nobasos_dynamicbins_", jmark, ".2022-04-12.rds"))
  
  LL.all <- readRDS(inf)
  
  LL.mat <- t(LL.all)
  
  # indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/repressive_cleaned"
  indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_cleaned_up_update_k27me3_umap"
  inf.meta <- file.path(indir.meta, paste0("metadata.", jmark, ".2022-03-01.txt"))
  dat.meta <- fread(inf.meta) %>%
    dplyr::select(-ctype.from.LL)

  # if (jmark %in% c("k4me1", "k4me3")){
  # } else if (jmark %in% c("k27me3", "k9me3")) {
  #   indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/repressive_cleaned"
  #   inf.meta <- file.path(indir.meta, paste0("metadata_celltyping_", jmark, ".dynamicbins.2022-02-18.txt"))
  #   dat.meta <- fread(inf.meta)
  # }
  ctypes.indx <- colnames(LL.mat)
  LL.mat.progs <- LL.mat[, prog.ctypes]
  ctypes.indx.progs <- colnames(LL.mat.progs)
  
  LL.best <- apply(LL.mat, 1, function(jrow){
    ctypes.indx[which.max(jrow)]
  })
  
  LL.secondbest <- apply(LL.mat, 1, function(jrow){
    names(jrow) <- colnames(LL.mat)
    jrow.nomax <- jrow[jrow!=max(jrow)]
    names(jrow.nomax)[which.max(jrow.nomax)]
  })
  
  LL.progs <- apply(LL.mat.progs, 1, function(jrow){
    ctypes.indx.progs[which.max(jrow)]
  })
  
  dat.LL.best <- data.frame(cell = names(LL.best), ctype.from.LL = LL.best, 
                            ctype.from.LL.secondbest = LL.secondbest, 
                            ctype.from.LL.prog = LL.progs, 
                            stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(ctype.from.LL = RemoveOldFromCtype(ctype.from.LL), 
           ctype.from.LL.secondbest = RemoveOldFromCtype(ctype.from.LL.secondbest), 
           ctype.from.LL.prog = RemoveOldFromCtype(ctype.from.LL.prog))
  
  dat.meta.reannotate <- dat.meta %>%
    left_join(., dat.LL.best)
  
  # if (jmark == "k4me3"){
  #   dat.meta.reannotate <- subset(dat.meta.reannotate, ! louvain %in% c("3", "7"))
  # }


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
  
  # dat.meta.reannotate <- dat.meta.reannotate %>%
  #     rowwise() %>%
  #   mutate(ctype.from.LL = ifelse(ctype.from.LL == "Basophils" & ctype != "Basophils", ctype.from.LL.secondbest, ctype.from.LL))
  # print(table(dat.meta.reannotate$ctype.from.LL))
  
  # if (jmark  "k4me3"){
    
    print("Adjs for all")
    # changing baso to ctype causes some to be annotated as HSPCs or DCs, which need to be changed
    # dat.meta.reannotate <- dat.meta.reannotate %>%
    #    rowwise() %>%
    #    mutate(ctype.from.LL.tmp = ifelse(ctype.from.LL == "Basophils" & ctype != "Basophils" & batch == "Old", ctype, NA),
    #           ctype.from.LL.tmp = gsub("^HSPCs$", "HSCs", ctype.from.LL.tmp),   
    #           ctype.from.LL.tmp = gsub("^DCs$", "Monocytes", ctype.from.LL.tmp), 
    #           ctype.from.LL = ifelse(!is.na(ctype.from.LL.tmp), ctype.from.LL.tmp, ctype.from.LL))
    # 
    # dat.meta.reannotate <- dat.meta.reannotate %>%
    #   rowwise() %>%
    #   mutate(ctype.from.LL = ifelse(ctype == "Basophils" & batch == "Old", ctype, ctype.from.LL))
    
    
  # } 
 
  # adjust pDCs  and NKs only in k27me3
  if (jmark == "k27me3"){
    dat.meta.reannotate <- dat.meta.reannotate %>%
        rowwise() %>%
        mutate(ctype.from.LL = ifelse(ctype.from.LL == "pDCs" & ctype != "pDCs", ctype.from.LL.secondbest, ctype.from.LL))
               # ctype.from.LL = ifelse(ctype.from.LL == "Basophils", ctype.from.LL.progs, ctype.from.LL))
  }
    

# End adjusts -------------------------------------------------------------
  
  # dat.meta.reannotate <- subset(dat.meta.reannotate, select = -ctype.from.LL.secondbest)
  print(table(dat.meta.reannotate$ctype.from.LL))
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
  
  m1.final.louv <- ggplot(dat.meta.reannotate, aes(x = umap1, y = umap2, color = colcode)) +
    geom_point() + 
    theme_bw() + 
    facet_wrap(~louvain) + 
    ggtitle(paste(jmark, "New + Old")) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m1.final.louv)
  
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







