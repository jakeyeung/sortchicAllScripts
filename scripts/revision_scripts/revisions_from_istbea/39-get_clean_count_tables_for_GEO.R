# Jake Yeung
# Date of Creation: 2022-05-24
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/39-get_clean_count_tables_for_GEO.R
# Processed data 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks


# Load metadata  ----------------------------------------------------------

inf.metas <- lapply(jmarks, function(jmark){
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/K562/count_tables/K562_", jmark, "/qc_metadata_new_only.K562_", jmark, ".0.8_0.5_3000.2022-03-11.txt")
  return(inf.meta)
})

dat.metas <- lapply(inf.metas, fread)


# Load count tables  ------------------------------------------------------

inf.mats <- lapply(jmarks, function(jmark){
  inf.mat <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/K562/count_tables/K562_", jmark, "/count_mat_new_only.K562_", jmark, ".0.8_0.5_3000.2022-03-11.rds")
  return(inf.mat)
})

dat.mats <- lapply(inf.mats, readRDS)

# Filter count tables -----------------------------------------------------

dat.mats.filt <- lapply(jmarks, function(jmark){
  jsub <- subset(dat.metas[[jmark]], is.good)
  cells.keep <- jsub$samp
  dat.mats[[jmark]][, cells.keep]
})


# Write outputs -----------------------------------------------------------



# Load bone marrows -------------------------------------------------------


countmat.bm <- lapply(jmarks, function(jmark){
  print(jmark)
  if (jmark == "k4me1"){
    inf.ldaout <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_varfilt/ldaAnalysis_fripfilt_varfilt/lda_outputs.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.k4me1.2022-01-28/ldaOut.count_mat_var_filt_dynamicbins.out_dynamic_bins_new_only.varcutoff.k4me1.2022-01-28.Robj"
  } else if (jmark == "k4me3"){
    inf.ldaout <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k4me3_cleaned/lda_outputs.count_mat_cleaned_no3no7_dynbins_allcells.k4me3.2022-04-12/ldaOut.count_mat_cleaned_no3no7_dynbins_allcells.k4me3.2022-04-12.Robj"
  } else if (jmark == "k27me3"){
    inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k27me3_clean_eryths/lda_outputs.count_mat_merged_with_old_dynbins.k27me3.2022-04-15/ldaOut.count_mat_merged_with_old_dynbins.", jmark, ".2022-04-15.Robj")
  } else if (jmark == "k9me3"){
    inf.ldaout <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_repressive_cleaned_from_jupyter/ldaAnalysis_fripfilt_varfilt_binfilt/lda_outputs.count_mat_cleaned_dynbins.", jmark, ".2022-02-16/ldaOut.count_mat_cleaned_dynbins.", jmark, ".2022-02-16.Robj")
  }
  load(inf.ldaout, v=T)
  return(count.mat)
})


indir.meta.bm.clean <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/ctypes_on_umap_batch_corrected_colors_fixed"

dat.metas.lsk.bm <- lapply(jmarks, function(jmark){
  inf.meta.bm <- file.path(indir.meta.bm.clean, paste0("LSK_metadata_colorcode_more_contrast.", jmark, ".2022-05-17.txt"))
  fread(inf.meta.bm)
})

dat.metas.umap.bm <- lapply(jmarks, function(jmark){
  inf.meta.bm <- file.path(indir.meta.bm.clean, paste0("umap_metadata_color_DC_monocyte_fixed.", jmark, ".2022-05-17.txt"))
  fread(inf.meta.bm)
})



# Save outputs ------------------------------------------------------------

# outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databank_GEO/processed"
outdir <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databank_GEO/processed_", Sys.Date())
dir.create(outdir)

for (jmark in jmarks){
  print(jmark)
  
  outf.meta.k562 <- file.path(outdir, paste0("metadata_K562_", jmark, ".txt"))
  outf.mat.k562 <- file.path(outdir, paste0("countmat_K562_", jmark, ".txt"))
  outf.meta.bm <- file.path(outdir, paste0("metadata_bonemarrow_allmerged_", jmark, ".txt"))
  outf.lsk.bm <- file.path(outdir, paste0("metadata_bonemarrow_LSKstained_", jmark, ".txt"))
  outf.mat.bm <- file.path(outdir, paste0("countmat_bonemarrow_allmerged_", jmark, ".txt"))
  
  fwrite(dat.metas[[jmark]], file = outf.meta.k562, sep = "\t")
  write.table(as.matrix(dat.mats[[jmark]]), file = outf.mat.k562, sep = "\t", row.names = TRUE, col.names = NA)
  
  fwrite(dat.metas.umap.bm[[jmark]], file = outf.meta.bm, sep = "\t")
  fwrite(dat.metas.lsk.bm[[jmark]], file = outf.lsk.bm, sep = "\t")
  write.table(as.matrix(countmat.bm[[jmark]]), file = outf.mat.bm, sep = "\t", row.names = TRUE, col.names = NA)
}

