# Jake Yeung
# Date of Creation: 2022-02-18
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/16-visualize_trajectories.R
# 


rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(topicmodels)

jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 30
jsettings[["min_dist"]] <- 0.1
jsettings[["random_state"]] <- 123


jmarks <- c("k4me1", "k4me3", "k27me3"); names(jmarks) <- jmarks

# jctypes <- c("Basophils", "Bcells", "DCs", "Eryths", "Granulocytes", "Monocytes", "NKs", "pDCs")
jctypes <- c("Bcells", "DCs", "Eryths", "Granulocytes", "Monocytes", "NKs", "pDCs")
names(jctypes) <- jctypes

# jctype <- "Bcells"
# jmark <- "k27me3"
jmark <- "k4me3"

jsuffix <- "oldandnew"
outdir <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs_plots/suffix_", jsuffix)
dir.create(outdir)
for (jmark in jmarks){
  
  
  print(jmark)
  
  
  
  for (jctype in jctypes){
    
    
    # jctype <- "DCs"
    # jmark <- "k4me1"
    
    # Load meta ---------------------------------------------------------------
    
    
    inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_cleaned_up/metadata.", jmark, ".txt")
    dat.meta <- fread(inf.meta)
    dat.meta.filt <- dat.meta %>%
      dplyr::select(cell, ctype.from.LL, colcode, batch)
    
    # 
    # # Load LDA ----------------------------------------------------------------
    # 
    # 
    # inf.lda <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_cluster/LDA_outputs_split_by_trajs/lda_trajs.countmat_tss.", jctype, ".", jmark, ".2022-02-17/ldaOut.countmat_tss.Monocytes.", jmark, ".2022-02-17.Robj")
    # load(inf.lda, v=T)
    # 
    # tm <- posterior(out.lda)
    # 
    # dat.umap <- DoUmapAndLouvain(tm$topics, jsettings)
    # 
    # dat.umap.annot <- left_join(dat.umap, dat.meta.filt)
    # 
    # ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = colcode)) + 
    #   geom_point() + 
    #   theme_bw() + 
    #   scale_color_identity() + 
    #   ggtitle(paste(jmark, jctype)) + 
    #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    
    
    
    # Load glmpca -------------------------------------------------------------
    
    inf.glmpca <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_cluster/glmpca_outputs_split_by_trajs/glmpca.", jctype, ".", jmark, ".bincutoff_0.binskeep_0.platename_batch.szname_none.niter_1000.RData")
    load(inf.glmpca, v=T)
    
    dat.glm <- data.frame(cell = rownames(glm.out$factors), glm.out$factors, stringsAsFactors = FALSE) %>%
      left_join(., dat.meta.filt)
    
    outf <- file.path(outdir, paste0("trajs_plots.", jctype, ".", jmark, ".", Sys.Date(), ".pdf"))
    pdf(file = outf, useDingbats = FALSE)
    
    # Print all dims ----------------------------------------------------------
    
    dims <- seq(1, ncol(glm.out$factors) - 1)
    for (d1 in dims){
      print(d1)
      d2 <- d1 + 1
      d1name <- paste0("dim", d1)
      d2name <- paste0("dim", d2)
      m.tmp1 <- ggplot(dat.glm, aes_string(x = d1name, y = d2name, color = "colcode")) + 
        geom_point() + 
        theme_bw() + 
        facet_wrap(~batch) + 
        scale_color_identity() + 
        ggtitle(paste(jmark, jctype)) + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      print(m.tmp1)
      
      m.tmp <- ggplot(dat.glm, aes_string(x = d1name, y = d2name, color = "colcode")) + 
        geom_point() + 
        theme_bw() + 
        scale_color_identity() + 
        ggtitle(paste(jmark, jctype)) + 
        facet_wrap(~ctype.from.LL) + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      print(m.tmp)
    }
    
    dev.off()
  }
  
  
  
}

