# Jake Yeung
# Date of Creation: 2022-02-17
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/2-check_LDA_output_dynamic_bins_repressive_cleaned.R
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
jsettings[["n_neighbors"]] <- 75
jsettings[["min_dist"]] <- 0.9
jsettings[["spread"]] <- 8
jsettings[["random_state"]] <- 123


# jmark <- "k4me1"

jsuffix1 <- "cleaned"
jsuffixbins <- "dynamicbins"
# jsuffixbins <- "allbins"

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata/repressive_cleaned"
dir.create(outdir)

jmarks <- c("k27me3", "k9me3"); names(jmarks) <- jmarks

# for (jmark in jmarks){
indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_cleaned_up"

dat.meta.var.lst <- parallel::mclapply(jmarks, function(jmark){
  print(jmark)
  # inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/metadata_celltyping_", jmark, "_dynamic_bins_merged_with_old.2022-01-28.txt")
  inf.meta <- file.path(indir.meta, paste0("metadata.", jmark, ".txt"))
  outf.meta <- file.path(outdir, paste0("metadata_plate_experi_batch.", jsuffix1, ".", jsuffixbins, ".", jmark, ".", Sys.Date(), ".txt"))
  outf.pdf <- file.path(outdir, paste0("LDA_outputs", jsuffix1, ".",  jsuffixbins, ".", jmark, ".", Sys.Date(), ".pdf"))
  indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_repressive_cleaned_from_jupyter/ldaAnalysis_fripfilt_varfilt_binfilt"
  dname <- paste0("count_mat_cleaned_dynbins.", jmark, ".2022-02-16")
  inf <- file.path(indir, paste0("lda_outputs.", dname), paste0("ldaOut.", dname, ".Robj"))
  assertthat::assert_that(file.exists(inf))
  # inf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt_again/ldaAnalysis_fripfilt_varfilt_binfilt/lda_outputs.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-02-13/ldaOut.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-02-13.Robj")
  
  dat.meta.ctypes <- fread(inf.meta)  %>%
    dplyr::select(cell, plate, experi, batch, mark, Batch, platerow, platecol, ctype) %>%
    rowwise() %>%
    mutate(ctype = gsub("BCells", "Bcells", ctype))
  
  pdf(outf.pdf, useDingbats = FALSE)
  
  load(inf, v=T)
  
  tm <- posterior(out.lda)
  
  dat.umap <- DoUmapAndLouvain(tm$topics, jsettings = jsettings) %>%
    left_join(., dat.meta.ctypes)
  
  # add var
  jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
  jchromos <- jchromos[!grepl("chrX|chrY", jchromos)]
  
  dat.impute.log <- t(log2(tm$topics %*% tm$terms))
  print(dat.impute.log[1:5, 1:5])
  dat.var <- CalculateVarAll(dat.impute.log, jchromos)
  dat.meta.new.var <- left_join(dat.umap, dat.var)
  
  
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
  
  
  m <- ggplot(dat.meta.new.var, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_viridis_c(direction = -1) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
  # save output metafiles
  fwrite(dat.meta.new.var, file = outf.meta)
  
  dev.off()
  return(dat.meta.new.var)
}, mc.cores = length(jmarks))
# }

print("Done mclapply")
outf.rds <- file.path(outdir, paste0("metadata_plate_experi_batch.", jsuffix1, ".", jsuffixbins, ".allmarks.", Sys.Date(), ".rds"))
saveRDS(dat.meta.var.lst, file = outf.rds)

