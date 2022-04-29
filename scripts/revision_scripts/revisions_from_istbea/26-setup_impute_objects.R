# Jake Yeung
# Date of Creation: 2022-04-24
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/26-setup_impute_objects.R
# 



rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)


jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks  # bins
jmarks.sub <- c("k4me1", "k4me3", "k27me3"); names(jmarks.sub) <- jmarks.sub  # TSS


# Load LDAs (TSS or bins) ---------------------------------------------------------------


dat.impute.lst <- lapply(jmarks, function(jmark){
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
  tm <- posterior(out.lda)
  dat.impute <- log2(t(tm$topics %*% tm$terms))
  return(dat.impute)
})


dat.impute.tss.lst <- lapply(jmarks.sub, function(jmark){
  if (jmark == "k4me1"){
    inf.ldaout.tss <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_TSS/ldaAnalysis_fripfilt_TSS/ldaOut.count_mat_TSS_combined.k4me1.2022-02-07.Robj"
  } else if (jmark == "k4me3"){
    inf.ldaout.tss <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k4me3_cleaned/lda_outputs.count_mat_cleaned_no3no7_tss_allcells.k4me3.2022-04-12/ldaOut.count_mat_cleaned_no3no7_tss_allcells.k4me3.2022-04-12.Robj"
  } else if (jmark == "k27me3"){
    # new only
    inf.ldaout.tss <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_k27me3_clean_eryths/lda_outputs.count_mat_new_only_TSS.k27me3.2022-04-15/ldaOut.count_mat_new_only_TSS.k27me3.2022-04-15.Robj")
  }
  load(inf.ldaout.tss, v=T)
  tm <- posterior(out.lda)
  dat.impute <- log2(t(tm$topics %*% tm$terms))
  return(dat.impute)
})


# Replace with k27me3  ----------------------------------------------------


inf.impute.k27me3.bc <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/mat_wide_k27me3_batch_corrected.2022-04-19.rds"
dat.impute.k27me3.bc <- readRDS(inf.impute.k27me3.bc)

inf.impute.k27me3.tss.bc <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/mat_wide_k27me3_merged_with_old_TSS_marker_genes_batch_corrected.2022-04-24.rds"
dat.impute.k27me3.tss.bc <- readRDS(inf.impute.k27me3.tss.bc)

dat.impute.lst$k27me3 <- dat.impute.k27me3.bc
dat.impute.tss.lst$k27me3 <- dat.impute.k27me3.tss.bc




# Handle metas ------------------------------------------------------------



# Load metas --------------------------------------------------------------

indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping_update_ctypes_from_LDA_k4me3_cleaned_k27me3_eryths2"

dat.meta.lst <- lapply(jmarks, function(jmark){
  inf.meta <- file.path(indir.meta, paste0("metadata_reannotate_from_LLmat_fix_ctypes_by_batch_dynamicbins.", jmark, ".txt"))
  fread(inf.meta)
})

# for K27me3 we have a new UMAP 
inf.meta.k27me3.bc <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/dat_umap_annot_batch_corrected.k27me3.2022-04-21.txt"
dat.meta.k27me3.bc <- fread(inf.meta.k27me3.bc)
dat.meta.lst$k27me3 <- dat.meta.k27me3.bc

dat.meta.colors <- subset(dat.meta.k27me3.bc, select = c(ctype.from.LL, colcode))
dat.meta.colors <- dat.meta.colors[!duplicated(dat.meta.colors), ]



# Save objects ------------------------------------------------------------

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/batch_effect_corrections/primetime_objects"

# saveRDS(dat.impute.lst, file.path(outdir, paste0("dat_impute_bins_lst.", Sys.Date(), ".rds")))
# saveRDS(dat.impute.tss.lst, file.path(outdir, paste0("dat_impute_tss_lst.", Sys.Date(), ".rds")))
saveRDS(dat.meta.colors, file.path(outdir, paste0("dat_meta_colors.", Sys.Date(), ".rds")))
saveRDS(dat.meta.lst, file.path(outdir, paste0("dat_meta_lst.", Sys.Date(), ".rds")))

for (jmark in jmarks){
  print(jmark)
  fwrite(dat.meta.lst[[jmark]], file = file.path(outdir, paste0("meta_data_", jmark, ".", Sys.Date(), ".txt")), sep = "\t")
  
  saveRDS(dat.impute.lst[[jmark]], file.path(outdir, paste0("dat_impute_bins_", jmark, ".", Sys.Date(), ".rds")))
  if (jmark != "k9me3"){
    saveRDS(dat.impute.tss.lst[[jmark]], file.path(outdir, paste0("dat_impute_tss_", jmark, ".", Sys.Date(), ".rds")))
  }
}
