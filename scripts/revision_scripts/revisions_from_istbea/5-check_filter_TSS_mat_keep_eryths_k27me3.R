# Jake Yeung
# Date of Creation: 2022-04-13
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/5-check_filter_TSS_mat_keep_eryths_k27me3.R
# 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# Load TSS ----------------------------------------------------------------

jmarksnew <- c("k27me3"); names(jmarksnew) <- jmarksnew

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables/filtered_count_tables_for_LDA_no_frip_filter_on_eryths/BM_k27me3"
dir.create(outdir)
indir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mat_TSS"
indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables/filtered_count_tables_for_LDA_no_frip_filter_on_eryths/BM_k27me3"

infs.tss <- lapply(jmarksnew, function(jmark){
  if (jmark != "k27me3"){
    inf <- file.path(indir, paste0("count_mat_TSS_combined.", jmark, ".2022-02-07.rds"))
  } else {
    inf <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mat_TSS_rep2rep3only/count_mat_TSS_rep2rep3only.k27me3.2022-02-08.rds"
  }
})

infs.meta <- lapply(jmarksnew, function(jmark){
  inf <- file.path(indir.meta, paste0("qc_metadata_new_only.", jmark, ".2022-04-12.txt"))
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

mats.tss.filt.newonly <- lapply(jmarksnew, function(jmark){
  print(jmark)
  jmat <- mats.tss[[jmark]]
  cells.keep <- (dat.metas[[jmark]] %>% filter(is.good & batch == "New"))$cell
  jmat.filt <- jmat[, cells.keep]
  print(dim(jmat))
  print(dim(jmat.filt))
  return(jmat.filt)
})



for (jmark in jmarksnew){
  print(jmark)
  outf.tmp <- file.path(outdir, paste0("count_mat_merged_with_old_TSS.", jmark, ".2022-04-12.rds"))
  saveRDS(mats.tss.filt[[jmark]], file = outf.tmp)
}


