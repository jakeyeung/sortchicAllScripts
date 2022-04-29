# Jake Yeung
# Date of Creation: 2022-04-12
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/8c-H3K4me3_cleaned_LDA_reannotate.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)

library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 30
jsettings[["min_dist"]] <- 0.1
jsettings[["random_state"]] <- 123
jsettings[["spread"]] <- 8



# Downstream --------------------------------------------------------------

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata_k4me3_cleaned_up"

jsuffixs <- c("no3", "no7", "no3no7"); names(jsuffixs) <- jsuffixs

inf.meta.all <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_oldannot/metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.k4me3.2022-04-12.txt"
dat.meta.all <- fread(inf.meta.all)
dat.meta.colors <- subset(dat.meta.all, select = c(ctype.from.LL, colcode))
colors.hash <- hash::hash(dat.meta.colors$ctype.from.LL, dat.meta.colors$colcode)

jmark <- "k4me3"


# jsuffix <- "no3no7"

inmain <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k4me3_cleaned_marcin/ldaAnalysis_cleaned"
# for (jsuffix in jsuffixs){
parallel::mclapply(jsuffixs, function(jsuffix){
  
  outpdf <- file.path(outdir, paste0("umap_k4me3_louvain_filtered.", jsuffix, ".", Sys.Date(), ".pdf"))
  pdf(outpdf, useDingbats = FALSE)
  
  outmeta <- file.path(outdir, paste0("metadata_k4me3_louvain_filtered.", jsuffix, ".", Sys.Date(), ".txt"))
  
  inf <- file.path(inmain, paste0("lda_outputs.count_mat_cleaned_", jsuffix, ".k4me3.2022-04-11"), paste0("ldaOut.count_mat_cleaned_", jsuffix, ".k4me3.2022-04-11.Robj"))
  load(inf, v=T)
  
  tm <- posterior(out.lda) %>%
    AddTopicToTmResult()
  
  dat.umap <- DoUmapAndLouvain(topics.mat = tm$topics, jsettings = jsettings)
  
  
  # Add annotations ---------------------------------------------------------
  
  dat.umap.annot <- left_join(dat.umap, subset(dat.meta.all, select = c(cell, ctype, ctype.from.LL, batch))) %>%
    rowwise() %>%
    mutate(colcode = colors.hash[[ctype.from.LL]])
  
  m <- ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = colcode)) + 
    geom_point() + 
    ggtitle(paste("k4me3 louvain filt", jsuffix)) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode,
                          guide = "legend") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  

  # Show each celltype individually  ----------------------------------------
  
  jctypes <- sort(unique(dat.umap.annot$ctype.from.LL)); names(jctypes) <- jctypes
  jbatchs <- c("Old", "New"); names(jbatchs) <- jbatchs
  
  for (jctype in jctypes){
    for (jbatch in jbatchs){
      jsub <- dat.umap.annot %>% 
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
      JFuncs::multiplot(m1.new2, m1.new2.LL, cols = 2)
    }
  }
  
  dev.off()
  fwrite(dat.umap.annot, file = outmeta, sep = "\t")
}, mc.cores = length(jsuffixs))

