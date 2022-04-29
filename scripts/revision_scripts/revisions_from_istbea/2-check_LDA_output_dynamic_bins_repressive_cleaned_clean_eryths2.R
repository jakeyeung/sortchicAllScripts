# Jake Yeung
# Date of Creation: 2022-04-16
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/2-check_LDA_output_dynamic_bins_repressive_cleaned_clean_eryths2.R
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
jsettings[["n_neighbors"]] <- 60
jsettings[["min_dist"]] <- 0.5
jsettings[["spread"]] <- 4
jsettings[["random_state"]] <- 123


jsuffixbins <- "dynamicbins"

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata/repressive_cleaned_clean_eryths2"
dir.create(outdir)

jmarks <- c("k27me3"); names(jmarks) <- jmarks

# jmark <- "k27me3"
jsuffix2.vec <- c("new_only_dynbins", "merged_with_old_dynbins")

# for (jmark in jmarks){
indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_cleaned_up"

for (jmark in jmarks){
  dat.meta.var.lst <- parallel::mclapply(jsuffix2.vec, function(jsuffix2){
    print(jmark)
    inf.meta <- file.path(indir.meta, paste0("metadata.", jmark, ".txt"))
    outf.meta <- file.path(outdir, paste0("metadata_plate_experi_batch.", jsuffix2, ".", jsuffixbins, ".", jmark, ".", Sys.Date(), ".txt"))
    outf.pdf <- file.path(outdir, paste0("LDA_outputs_adj_uma", jsuffix2, ".",  jsuffixbins, ".", jmark, ".", Sys.Date(), ".pdf"))
    indir <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_", jmark, "_clean_eryths")
    dname <- paste0("count_mat_", jsuffix2, ".", jmark, ".2022-04-15")
    inf <- file.path(indir, paste0("lda_outputs.", dname), paste0("ldaOut.", dname, ".Robj"))
    assertthat::assert_that(file.exists(inf))
    
    dat.meta.ctypes <- fread(inf.meta)  %>%
      dplyr::select(cell, plate, experi, batch, mark, Batch, platerow, platecol, ctype) %>%
      rowwise() %>%
      mutate(ctype = gsub("BCells", "Bcells", ctype))
    
    pdf(outf.pdf, useDingbats = FALSE)
    
    load(inf, v=T)
    
    tm <- posterior(out.lda)
    
    dat.umap <- DoUmapAndLouvain(tm$topics, jsettings = jsettings) %>%
      left_join(dat.meta.ctypes) %>%
      LabelMetadataByCtype()
    
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
    
    
    m <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = experi)) + 
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
    
    m <- ggplot(dat.meta.new.var, aes(x = umap1, y = umap2, color = ctype)) + 
      geom_point() + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
    
    m <- ggplot(dat.meta.new.var, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
      geom_point() + 
      theme_bw() + 
      facet_wrap(~ctype) + 
      scale_color_viridis_c(direction = -1) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
    
    m <- ggplot(dat.meta.new.var, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
      geom_point() + 
      theme_bw() + 
      facet_wrap(~batch) + 
      scale_color_viridis_c(direction = -1) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
    
    
    # save output metafiles
    fwrite(dat.meta.new.var, file = outf.meta)
    
    dev.off()
    
    outf.rds <- file.path(outdir, paste0("metadata_plate_experi_batch.", jsuffix2, ".", jsuffixbins, ".allmarks.", Sys.Date(), ".rds"))
    saveRDS(dat.meta.new.var, file = outf.rds)
    
    return(dat.meta.new.var)
  }, mc.cores = length(jsuffix2.vec))
  # }

}

