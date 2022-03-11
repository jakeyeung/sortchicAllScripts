# Jake Yeung
# Date of Creation: 2022-02-14
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/5-check_var_filt_k27me3_only.R
# 

rm(list=ls()) 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("k4me1", "k4me3", "k9me3")
# jmark <- "k27me3"
outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mats_k27me3_var_filt"
dir.create(outdir)

for (jmark in jmarks){
  print(jmark)
  
  outpdf <- file.path(outdir, paste0("plots_varfilt_", jmark, ".", Sys.Date(), ".pdf"))
  
  pdf(file = outpdf, useDingbats = FALSE)
  
  # Load UMAP ---------------------------------------------------------------
  
  inf.prime <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_with_colors.", jmark, ".2022-02-13.txt")
  dat.meta.prime <- fread(inf.prime)
  
  
  
  # Add var -----------------------------------------------------------------
  
  # # only new
  # inf.meta.var <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mats_varfilt2/metadata_with_var.out_dynamic_bins_new_only.varcutoff.k27me3.2022-02-14.txt"
  # old and new, but old is duplicated
  inf.meta.var <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mats_varfilt_again/metadata_with_var.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-02-13.txt")
  dat.meta.var <- fread(inf.meta.var)
  
  print(length(unique(dat.meta.var$cell)))
  
  
  
  # Show umap with var  -----------------------------------------------------
  
  
  dat.meta.prime.var <- left_join(dat.meta.prime, subset(dat.meta.var, select = c(cell, is.good, cell.var.within.sum.norm)))
  
  m <- ggplot(dat.meta.prime.var, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_viridis_c(direction = -1) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.meta.prime.var, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + 
    theme_bw() + 
    facet_wrap(~batch) + 
    scale_color_viridis_c(direction = -1) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.meta.prime.var, aes(x = umap1, y = umap2, color =  is.good)) + 
    geom_point() + 
    theme_bw() + 
    # scale_color_viridis_c(direction = -1) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.meta.prime.var, aes(x = umap1, y = umap2, color =  is.good)) + 
    geom_point() + 
    theme_bw() + 
    facet_wrap(~batch) + 
    # scale_color_viridis_c(direction = -1) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  
  # Get cells keep ----------------------------------------------------------
  
  # keep all old cells, remove new that have low var
  
  cells.keep.old <- subset(dat.meta.prime.var, batch == "Old")$cell
  cells.keep.new <- subset(dat.meta.prime.var, batch == "New" & is.good)$cell
  cells.keep <- c(cells.keep.old, cells.keep.new)
  
  assertthat::assert_that(length(cells.keep) == length(unique(cells.keep)))
  
  # Load counts -------------------------------------------------------------
  
  
  # load dynbins
  # load all bins
  
  inf.dynbins <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables/filtered_count_tables_for_LDA_bugfixed_dynamic_bins/BM_", jmark, "/count_mat_merged_with_old.", jmark, ".2022-01-26.rds")
  inf.allbins <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables/filtered_count_tables_for_LDA_bugfixed/BM_", jmark, "/count_mat_merged_with_old.", jmark, ".2022-01-26.rds")
  assertthat::assert_that(file.exists(inf.dynbins))
  assertthat::assert_that(file.exists(inf.allbins))
  
  count.dynbins <- readRDS(inf.dynbins)
  count.allbins <- readRDS(inf.allbins)
  
  # filter bins
  jchromos.check <- sort(unique(sapply(rownames(count.allbins), function(x) strsplit(x, ":")[[1]][[1]])))
  print(jchromos.check)
  rnames.orig <- rownames(count.allbins)
  jchromos.auto <- paste("chr", c(seq(19), "X", "Y"), sep = "")
  jchromos.grep <- paste(paste0("^", jchromos.auto), collapse = "|")
  rnames.i <- grepl(jchromos.grep, rnames.orig)
  rnames.filt <- rnames.orig[rnames.i]
  print(dim(count.allbins))
  count.allbins <- count.allbins[rnames.filt, ]
  print(dim(count.allbins))
  
  # load TSS
  if (jmark == "k27me3"){
    inf.tss <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mat_TSS_rep2rep3only/count_mat_TSS_rep2rep3only.", jmark, ".2022-02-08.rds")
  } else {
    inf.tss <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mat_TSS/count_mat_TSS_combined.", jmark, ".2022-02-07.rds")
  }
  count.tss <- readRDS(inf.tss)
  
  
  
  # Filter counts  ----------------------------------------------------------
  
  
  count.mat.filt.dynbins <- count.dynbins[, cells.keep]
  count.mat.filt.allbins <- count.allbins[, cells.keep]
  count.mat.filt.tss <- count.tss[, cells.keep]
  
  
  # Save outputs ------------------------------------------------------------
  
  
  print("Dim before:")
  print(dim(count.dynbins))
  print(dim(count.allbins))
  print(dim(count.tss))
  print("Dim after dynbins:")
  print(dim(count.mat.filt.dynbins))
  print("Dim after allbins:")
  print(dim(count.mat.filt.allbins))
  print("Dim after tss:")
  print(dim(count.mat.filt.tss))
  
  # write varmeta
  fwrite(dat.meta.prime.var, file = file.path(outdir, paste0("metadata_with_varfilt2.", jmark, ".txt")), sep = "\t")
  # write filtered count tables
  saveRDS(count.mat.filt.dynbins, file = file.path(outdir, paste0("count_mat_varfilt2_dynamicbins.", jmark, ".rds")))
  saveRDS(count.mat.filt.allbins, file = file.path(outdir, paste0("count_mat_varfilt2_allbins.", jmark, ".rds")))
  saveRDS(count.mat.filt.allbins, file = file.path(outdir, paste0("count_mat_varfilt2_TSS.", jmark, ".rds")))
  
  dev.off()
  
  
  
  
}


print('Done')