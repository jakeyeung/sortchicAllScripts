# Jake Yeung
# Date of Creation: 2022-02-14
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/5-check_var_filt_both_mats_unique_cellnames_binfilt.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(topicmodels)


# jmark <- "k9me3"
indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/metadata"
outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mats_varfilt4"
dir.create(outdir)

jsuffix1 <- "_dynamic_bins"
jsuffix2 <- "_new_only"

jmarksnew <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarksnew) <- jmarksnew
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksold) <- jmarksnew
# varcutoff <- 0.4
varcutoff.lst <- c(1.5, 1.5, 2, 1); names(varcutoff.lst) <- jmarksnew
# varcutoff.lst <- c(1.5, 1.5, 1.9, 1); names(varcutoff.lst) <- jmarksnew
# varcutoff.lst <- c(0.4, 1, 0.5, 0.4); names(varcutoff.lst) <- jmarksnew
# varcutoff.lst <- c(0.6, 1.2, 0.9, 0.6); names(varcutoff.lst) <- jmarksnew

for (jmark in jmarksnew){
  
  varcutoff <- varcutoff.lst[[jmark]]
  
  print(jmark)
  
  jsuffixout <- paste0("out", jsuffix1, jsuffix2, ".varcutoff.", jmark, ".", Sys.Date())
  
  pdf(file.path(outdir, paste0("plots_with_var.", jsuffixout, ".", Sys.Date(), ".pdf")), useDingbats = FALSE)
  
  jmarkold <- jmarksold[[jmark]]
  
  # indir <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_bugfixed", jsuffix1)
  # inf <- file.path(indir, paste0("ldaAnalysis_fripfilt_BM_", jmark, "/lda_outputs.count_mat", jsuffix2, ".", jmark, ".2022-01-26/ldaOut.count_mat", jsuffix2, ".", jmark, ".2022-01-26.Robj"))
  inf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_allbins.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28/ldaOut.count_mat_var_filt_allbins.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28.Robj")
  assertthat::assert_that(file.exists(inf))
  load(inf, v=T)
  
  tm <- posterior(out.lda)
  
  # inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/metadata_celltyping_", jmark, "_dynamic_bins_merged_with_old.2022-01-28.txt")
  # inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/celltyping/metadata_celltyping_", jmark, "_dynamic_bins_merged_with_old.2022-01-28.txt")
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_with_colors.", jmark, ".2022-02-13.txt")
  dat.meta.new <- fread(inf.meta)
  dat.meta.old <- subset(dat.meta.new, batch == "Old")
  
  cells.from.old <- dat.meta.old$cell
  
  cells.from.mat <- colnames(count.mat)
  dat.meta.new <- subset(dat.meta.new, cell %in% cells.from.mat)  # all cells
  
  # load meta
  m <- ggplot(dat.meta.new, aes(x = umap1, y = umap2, color = ctype)) + 
    geom_point() + 
    facet_wrap(~experi) + 
    theme_bw() + 
    ggtitle(paste(jmark, jsuffix1, jsuffix2)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  m <- ggplot(dat.meta.new, aes(x = umap1, y = umap2, color = experi)) + 
    geom_point() + 
    facet_wrap(~ctype) + 
    theme_bw() + 
    ggtitle(paste(jmark, jsuffix1, jsuffix2)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
  print(m)
  
  # check var
  
  # jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
  dat.impute.log <- t(log2(tm$topics %*% tm$terms))
  # jchromos <- sort(unique(sapply(rownames(dat.impute.log), function(x) strsplit(x, ":")[[1]][[1]])))
  # jchromos <- jchromos[!grepl("chrX|chrY", jchromos)]
  jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
  jchromos <- jchromos[!grepl("chrX|chrY", jchromos)]
  print(dat.impute.log[1:5, 1:5])
  dat.var <- CalculateVarAll(dat.impute.log, jchromos)
  
  dat.meta.new.var <- left_join(dat.meta.new, dat.var)
  
  m <- ggplot(dat.meta.new.var, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_viridis_c(direction = -1) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  # cutoff
  dat.meta.new.var.cutoff <- dat.meta.new.var %>%
    rowwise() %>%
    mutate(is.good = cell.var.within.sum.norm >= varcutoff)
  
  m <- ggplot(dat.meta.new.var.cutoff, aes(x = umap1, y = umap2, color = is.good)) + 
    geom_point() + 
    ggtitle(jmark, paste("Cutoff:", varcutoff)) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.meta.new.var.cutoff %>% filter(is.good), aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_viridis_c(direction = -1) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.meta.new.var.cutoff %>% filter(batch == "Old" | is.good), aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle("Old batch, is good") + 
    scale_color_viridis_c(direction = -1) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(dat.meta.new.var.cutoff, aes(x = cell.var.within.sum.norm)) + 
    geom_density() + 
    geom_vline(xintercept = varcutoff) + 
    theme_bw() + 
    ggtitle(jmark, paste("Cutoff:", varcutoff)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  dev.off()
  
  # new batch only for kee pnew
  cells.keep.new <- subset(dat.meta.new.var.cutoff, is.good & batch == "New")$cell
  cells.keep <- c(cells.from.old, cells.keep.new)
  assertthat::assert_that(length(cells.keep) == length(unique(cells.keep)))
  
  # check we have both old and new
  experis.new <- sapply(cells.keep.new, function(x) ClipLast(x, jsep = "_"))
  print(unique(experis.new))
  
  experis.old <- sapply(cells.from.old, function(x) ClipLast(x, jsep = "_"))
  print(unique(experis.old))
  
  experis <- sapply(cells.keep, function(x) ClipLast(x, jsep = "_"))
  print(unique(experis))
  
  
  # load countmats: dynamic bins
  inf.dynbins <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables/filtered_count_tables_for_LDA_bugfixed_dynamic_bins/BM_", jmark, "/count_mat_merged_with_old.", jmark, ".2022-01-26.rds")
  inf.allbins <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables/filtered_count_tables_for_LDA_bugfixed/BM_", jmark, "/count_mat_merged_with_old.", jmark, ".2022-01-26.rds")
  assertthat::assert_that(file.exists(inf.dynbins))
  assertthat::assert_that(file.exists(inf.allbins))
  
  count.dynbins <- readRDS(inf.dynbins)
  count.allbins <- readRDS(inf.allbins)
  
  # filter allbins
  # check before
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
  
  jchromos.check <- sort(unique(sapply(rownames(count.allbins), function(x) strsplit(x, ":")[[1]][[1]])))
  print(jchromos.check)
  
  
  # load countmats: merge with old 
  
  count.mat.filt.dynbins <- count.dynbins[, cells.keep]
  count.mat.filt.allbins <- count.allbins[, cells.keep]
  
  print("Dim before:")
  print(dim(count.dynbins))
  print(dim(count.allbins))
  print("Dim after dynbins:")
  print(dim(count.mat.filt.dynbins))
  print("Dim after allbins:")
  print(dim(count.mat.filt.allbins))
  
  # write varmeta
  fwrite(dat.meta.new.var.cutoff, file = file.path(outdir, paste0("metadata_with_var.", jsuffixout, ".txt")), sep = "\t")
  # write filtered count tables
  saveRDS(count.mat.filt.dynbins, file = file.path(outdir, paste0("count_mat_var_filt_dynamicbins.", jsuffixout, ".rds")))
  saveRDS(count.mat.filt.allbins, file = file.path(outdir, paste0("count_mat_var_filt_allbins.", jsuffixout, ".rds")))
}

  
  
  
  
  
 

