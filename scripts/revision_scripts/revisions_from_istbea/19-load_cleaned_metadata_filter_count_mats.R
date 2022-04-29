# Jake Yeung
# Date of Creation: 2022-03-22
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/19-load_cleaned_metadata_filter_count_mats.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


jmarksnew <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarksnew) <- jmarksnew
jmarksrepress <- c("k27me3", "k9me3"); names(jmarksrepress) <- jmarksrepress

# Load metas for cell types --------------------------------------------------------------

indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/metadata_cleaned_up"

dat.meta.lst <- lapply(jmarksnew, function(jmark){
  fread(file.path(indir.meta, paste0("metadata.", jmark, ".txt")))
})


# Update UMAPs with re-LDA'd repressive marks -----------------------------

indir.cleaned.meta <- "metadata_plate_experi_batch.cleaned.dynamicbins.k27me3.2022-02-17.txt"

dat.cleaned.meta.lst <- lapply(jmarksnew, function(jmark){
  if (!jmark %in% jmarksrepress){
    dat.meta <- dat.meta.lst[[jmark]]
  } else {
    inf.meta.cleaned <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/metadata/repressive_cleaned/metadata_plate_experi_batch.cleaned.dynamicbins.", jmark, ".2022-02-17.txt")
    dat.meta.cleaned <- fread(inf.meta.cleaned) %>%
      dplyr::select(cell, umap1, umap2) 
    assertthat::assert_that(nrow(dat.meta.cleaned) > 0)
    dat.meta.old <- dat.meta.lst[[jmark]] %>%
      dplyr::select(-umap1, -umap2)
    dat.meta <- left_join(dat.meta.old, dat.meta.cleaned)
  }
  return(dat.meta)
})

jalpha <- 0.5
m.lst <- lapply(jmarksnew, function(jmark){
  jdat <- dat.cleaned.meta.lst[[jmark]]
  ncells <- nrow(jdat)
  jsub <- jdat %>%
    dplyr::select(ctype.from.LL, colcode) 
  jsub <- jsub[!duplicated(jsub), ]
  
  ggplot(jdat, aes(x = umap1, y = umap2, color = colcode)) + 
    geom_point(alpha = jalpha) + 
    theme_bw() + 
    ggtitle(jmark, paste0("ncells: ", ncells)) + 
    scale_color_identity( labels = jsub$ctype.from.LL, breaks = jsub$colcode,
                          guide = "legend")  + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

print(m.lst)


# Load count data  --------------------------------------------------------

# indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mats_k27me3_var_filt"
indir.bins <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mats_varfilt"
indir.tss <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mat_TSS"

count.mats.allbins <- lapply(jmarksnew, function(jmark){
  # inf.tmp <- file.path(indir, paste0("count_mat_varfilt2_allbins.", jmark, ".rds"))
  inf.tmp <- file.path(indir.bins, paste0("count_mat_var_filt_allbins.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28.rds"))
  readRDS(inf.tmp)
})

count.mats.dynamicbins <- lapply(jmarksnew, function(jmark){
  # inf.tmp <- file.path(indir, paste0("count_mat_varfilt2_dynamicbins.", jmark, ".rds"))
  inf.tmp <- file.path(indir.bins, paste0("count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-01-28.rds"))
  readRDS(inf.tmp)
})

count.mats.tss <- lapply(jmarksnew, function(jmark){
  # inf.tmp <- file.path(indir, paste0("count_mat_varfilt2_TSS.", jmark, ".rds"))
  inf.tmp <- file.path(indir.tss, paste0("count_mat_TSS_combined.", jmark, ".2022-02-07.rds"))
  readRDS(inf.tmp)
})

# filter out 
count.mats.allbins.filt <- lapply(jmarksnew, function(jmark){
  cellskeep <- subset(dat.cleaned.meta.lst[[jmark]], batch == "New")$cell
  jmat <- count.mats.allbins[[jmark]][, cellskeep]
  return(jmat)
})

count.mats.dynamicbins.filt <- lapply(jmarksnew, function(jmark){
  cellskeep <- subset(dat.cleaned.meta.lst[[jmark]], batch == "New")$cell
  jmat <- count.mats.dynamicbins[[jmark]][, cellskeep]
  return(jmat)
})

count.mats.tss.filt <- lapply(jmarksnew, function(jmark){
  cellskeep <- subset(dat.cleaned.meta.lst[[jmark]], batch == "New")$cell
  jmat <- count.mats.tss[[jmark]][, cellskeep]
  return(jmat)
})



# Filter by new only ------------------------------------------------------

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mat_newonly_cleaned"

for (jmark in jmarksnew){
  print(jmark)
  outbins <- file.path(outdir, paste0("countmat_cleaned_newonly_allbins.", jmark, ".", Sys.Date(), ".rds"))
  outdynbins <- file.path(outdir, paste0("countmat_cleaned_newonly_dynbins.", jmark, ".", Sys.Date(), ".rds"))
  outtss <- file.path(outdir, paste0("countmat_cleaned_newonly_tss.", jmark, ".", Sys.Date(), ".rds"))
  
  outmeta <- file.path(outdir, paste0("metadata_cleaned_newonly.", jmark, ".", Sys.Date(), ".txt"))
  
  dat.meta.tmp <- subset(dat.cleaned.meta.lst[[jmark]], batch == "New")
  
  matbins <- count.mats.allbins.filt[[jmark]]
  matdynbins <- count.mats.dynamicbins.filt[[jmark]]
  mattss <- count.mats.tss.filt[[jmark]]
 
  print(dat.meta.tmp) 
  print(dim(matbins))
  print(dim(matdynbins))
  print(dim(mattss))
  
  saveRDS(matbins, file = outbins)
  saveRDS(matdynbins, file = outdynbins)
  saveRDS(mattss, file = outtss)
  
  fwrite(dat.meta.tmp, outmeta, sep = "\t")
}




