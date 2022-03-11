# Jake Yeung
# Date of Creation: 2022-02-14
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/2-check_LDA_output_dynamic_bins_varfilt2.R
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
library(scchicFuncs)


jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 50
jsettings[["min_dist"]] <- 0.9
jsettings[["spread"]] <- 8
jsettings[["random_state"]] <- 123


jmark <- "k4me1"

jsuffix1 <- "varfilt"
jsuffixbins <- "dynamicbins"
# jsuffixbins <- "allbins"

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata/varfilt_again"
dir.create(outdir)

jmarks <- c("k4me3", "k27me3", "k9me3")

for (jmark in jmarks){
  print(jmark)
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/metadata_celltyping_", jmark, "_dynamic_bins_merged_with_old.2022-01-28.txt")
  outf.meta <- file.path(outdir, paste0("metadata_plate_experi_batch.", jsuffix1, ".", jsuffixbins, ".", jmark, ".", Sys.Date(), ".txt"))
  outf.pdf <- file.path(outdir, paste0("LDA_outputs", jsuffix1, ".",  jsuffixbins, ".", jmark, ".", Sys.Date(), ".pdf"))
  # indir <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_", jsuffix1)
  # dname <- paste0("count_mat_var_filt_", jsuffixbins, ".out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28")
  # inf <- file.path(indir, "ldaAnalysis_fripfilt_varfilt", paste0("lda_outputs.", dname), paste0("ldaOut.", dname, ".Robj"))
  inf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt_again/ldaAnalysis_fripfilt_varfilt_binfilt/lda_outputs.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-02-13/ldaOut.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-02-13.Robj")
  assertthat::assert_that(file.exists(inf))
  
  dat.meta.ctypes <- fread(inf.meta)  %>%
    dplyr::select(cell, plate, experi, batch, mark, Batch, platerow, platecol, ctype) %>%
    rowwise() %>%
    mutate(ctype = gsub("BCells", "Bcells", ctype))
  
  pdf(outf.pdf, useDingbats = FALSE)
  
  load(inf, v=T)
  
  tm <- posterior(out.lda)
  
  dat.umap <- DoUmapAndLouvain(tm$topics, jsettings = jsettings) %>%
    left_join(., dat.meta.ctypes)
  
  m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
  m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = plate)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
  m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jmark) + 
    facet_wrap(~experi) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
  
  
  # save output metafiles
  fwrite(dat.umap, file = outf.meta)
  
  dev.off()
  
  
  
  
}
