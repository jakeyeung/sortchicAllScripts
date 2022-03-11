t# Jake Yeung
# Date of Creation: 2022-02-13
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/5-check_filter_TSS_mat_by_varfilt_cells.R
# 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# Load TSS ----------------------------------------------------------------

jmarksnew <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarksnew) <- jmarksnew

indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mat_TSS"
indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mats_varfilt_again"

infs.tss <- lapply(jmarksnew, function(jmark){
  if (jmark != "k27me3"){
    inf <- file.path(indir, paste0("count_mat_TSS_combined.", jmark, ".2022-02-07.rds"))
  } else {
    inf <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mat_TSS_rep2rep3only/count_mat_TSS_rep2rep3only.k27me3.2022-02-08.rds"
  }
})

infs.meta <- lapply(jmarksnew, function(jmark){
  inf <- file.path(indir.meta, paste0("metadata_with_var.out_dynamic_bins_new_only.varcutoff.", jmark, ".2022-02-13.txt"))
})

mats.tss <- lapply(infs.tss, function(jinf){
  readRDS(jinf)
})

dat.metas <- lapply(infs.meta, function(jinf){
  fread(jinf)
})

# Filter by cells  --------------------------------------------------------

mats.tss.filt <- lapply(jmarksnew, function(jmark){
  print(jmark)
  jmat <- mats.tss[[jmark]]
  cells.keep <- (dat.metas[[jmark]] %>% filter(is.good))$cell
  jmat.filt <- jmat[, cells.keep]
  print(dim(jmat))
  print(dim(jmat.filt))
  return(jmat.filt)
})

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mat_TSS_varfilt_again"

for (jmark in jmarks){
  print(jmark)
  outf.tmp <- file.path(outdir, paste0("count_mat_TSS_combined_varfilt_again.", jmark, ".", Sys.Date(), ".rds"))
  saveRDS(mats.tss.filt[[jmark]], file = outf.tmp)
}


