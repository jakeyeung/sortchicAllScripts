# Jake Yeung
# Date of Creation: 2022-04-25
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/11-pseudotime_downstream_k27me3_k9me3_together_cleaned.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)

library(topicmodels)
library(DescTools)

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/ptimefits_hits_k27_k9_together_cleaned"
dir.create(outdir)


jmarks <- c("k27me3", "k9me3"); names(jmarks) <- jmarks

# load metadata -----------------------------------------------------------

jctype <- "Granulocytes"

jctypes <- c("Eryths", "Bcells", "Granulocytes"); names(jctypes) <- jctypes

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
for (jctype in jctypes){
  
  # inf.traj.k27 <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/by_ctype/trajs_outputs.k27me3.2022-02-20.rds")
  inf.traj.k27 <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/cleaned2_batch_corrected_eryth_fix/trajs_outputs_batch_corrected.k27me3.2022-04-21.rds")
  dat.traj.k27 <- readRDS(inf.traj.k27)
  
  inf.traj.k9 <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/cleaned/trajs_outputs.k9me3.rds")
  dat.traj.k9 <- readRDS(inf.traj.k9)
  
  dat.traj.lst <- list(dat.traj.k27[[jctype]] %>% filter(!is.na(ptime)), dat.traj.k9[[jctype]] %>% filter(!is.na(ptime)))
  names(dat.traj.lst) <- jmarks
  
  
  # Load imputes ------------------------------------------------------------
  
  # mat.impute.log.lst <- lapply(jmarks, function(jmark){
  #   inf.lda <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_allbins.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28/ldaOut.count_mat_var_filt_allbins.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28.Robj")
  #   load(inf.lda, v=T)
  #   tm <- posterior(out.lda)
  #   mat.impute.log <- t(log2(tm$topics %*% tm$terms))
  # })
  
  mat.impute.log.lst <- lapply(jmarks, function(jmark){
    inf.impute <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/primetime_objects/dat_impute_bins_", jmark, ".2022-04-24.rds")
    readRDS(inf.impute)
  })
  
  
  # Load metas --------------------------------------------------------------
  
  
  # indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_cleaned_up_update_k27me3_umap"
  indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/primetime_objects"
  dat.meta.full.lst <- lapply(jmarks, function(jmark){
    fread(file.path(indir.meta, paste0("meta_data_", jmark, ".2022-04-24.txt")))
  })
  dat.meta.lst <- lapply(jmarks, function(jmark){
    fread(file.path(indir.meta, paste0("meta_data_", jmark, ".2022-04-24.txt"))) %>%
      dplyr::select(cell, colcode, batch)
  })
  
  # ctype.from.LL to colcode hash
  ctype2col.hash <- hash::hash(dat.meta.full.lst$k27me3$ctype.from.LL, dat.meta.full.lst$k27me3$colcode)
  
  
  # Plot some hits ----------------------------------------------------------
  
  inf.hits.k9 <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/ptimefits_plots_top_hits/fits_with_annotations.", jctype, ".k9me3.2022-02-06.neg_slope.txt")
  inf.hits.k27 <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/ptimefits_plots_top_hits/fits_with_annotations.", jctype, ".k27me3.2022-02-05.neg_slope.txt")
  
  dat.hits.k9 <- fread(inf.hits.k9) %>%
    mutate(mark = "k9me3")
  dat.hits.k27 <- fread(inf.hits.k27) %>%
    mutate(mark = "k27me3")
  
  dat.hits.merge <- rbind(dat.hits.k9, dat.hits.k27)
  
  jhit <- "chr11:86800000-86850000"
  
  jhits.k27 <- dat.hits.k27$region_coord[1:15]
  jhits.k9 <- dat.hits.k9$region_coord[1:15]
  
  m.dens <- ggplot(dat.hits.merge, aes(x = estimate, fill = mark)) + 
    geom_density(alpha = 0.5) + 
    ggtitle(jctype) +  
    geom_vline(xintercept = 0, linetype = "dotted") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  print(m.dens)
  
  
  jhits.negslope <- c(jhits.k27, jhits.k9)
  
  outpdf <- file.path(outdir, paste0("traj_batchcorrect_by_ctype_cleaned_", jctype, ".neg_slope_hits.", Sys.Date(), ".pdf"))
  
  pdf(outpdf, useDingbats = FALSE)
  print(m.dens)
  for (jhit in jhits.negslope){
    
    jgene <- unique(subset(dat.hits.merge, region_coord == jhit)$gene)
    print(jgene)
    
    if (length(jgene) > 0){
      # choose one
      jgene <- jgene[[1]]
    }
    
    jtitle <- paste(jgene, jhit, sep = ",")
    print(jtitle)
    
    dat.exprs.long <- lapply(jmarks, function(jmark){
      jmat <- mat.impute.log.lst[[jmark]]
      dat.long <- data.frame(signal = jmat[jhit, ], cell = colnames(jmat), stringsAsFactors = FALSE) %>%
        left_join(dat.traj.lst[[jmark]], .) %>%
        rowwise() %>%
        mutate(colcode = ctype2col.hash[[ctype.from.LL]]) %>%
        # left_join(., dat.meta.lst[[jmark]]) %>%
        mutate(mark = jmark) %>%
        group_by(mark, batch) %>%
        # filter(batch == "New") %>%
        mutate(ptime.norm = (ptime - min(ptime)) / (max(ptime) - min(ptime)))
      
      # # batch corret
      # dat.long <- dat.long %>%
      #   group_by(mark) %>%
      #   mutate(ptime.norm = )
      
    }) %>%
      bind_rows()
    
    
    # cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
    m0 <- ggplot(dat.exprs.long, aes(x = ptime.norm, y = signal, color = mark)) + 
      geom_point(alpha = 0.5) + 
      ggtitle(jtitle) + 
      theme_bw() + 
      scale_color_manual(values = cbPalette) + 
      theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    # 
    m1 <- ggplot(dat.exprs.long, aes(x = ptime.norm, y = signal, color = colcode, shape = mark)) + 
      geom_point() + 
      theme_bw() + 
      ggtitle(jtitle) + 
      scale_color_identity() + 
      theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    # 
    # m2 <- ggplot(dat.exprs.long, aes(x = ptime.norm, y = signal, color = colcode, shape = mark)) + 
    #   geom_point() + 
    #   theme_bw() + 
    #   ggtitle(jtitle) + 
    #   facet_wrap(~mark, ncol = 1) + 
    #   scale_color_identity() + 
    #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    # 
    print(m0)
    print(m1)
    # print(m2)
  }
  dev.off()
  
  
  
}


