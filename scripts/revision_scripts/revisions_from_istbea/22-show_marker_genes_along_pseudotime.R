# Jake Yeung
# Date of Creation: 2022-04-17
# File: ~/projects/scchic/scripts/revision_scripts/revisions_from_istbea/22-show_marker_genes_along_pseudotime.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)


jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks  # bins
jmarks.sub <- c("k4me1", "k4me3", "k27me3"); names(jmarks) <- jmarks  # TSS


# Load LDAs (TSS or bins) ---------------------------------------------------------------


tm.lst <- lapply(jmarks, function(jmark){
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
  return(tm)
})


tm.tss.lst <- lapply(jmarks.sub, function(jmark){
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
  return(tm)
})




# Load trajs --------------------------------------------------------------

indir.traj <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/trajs/cleaned"
dat.trajs <- lapply(jmarks, function(jmark){
  inf.trajs <- file.path(indir.traj, paste0("trajs_outputs.", jmark, ".rds"))
  readRDS(inf.trajs)
})



# Plot TSS  ---------------------------------------------------------------

markergenes <- c("Ly6c2", "Ly6g", "S100a8", "S100a2", "Chil3", "Sox6", "Tal1", "Gata1", "Ebf1", "Cd79a", "Cd79b", "Hoxa9", "Meis1", "Runx2", "Kit", "Hlf", "Erg", "Cd34", "Stat4", "Tcf7", "Cebpe")

jctypes <- c("Granulocytes", "Eryths", "Bcells"); names(jctypes) <- jctypes

# inf.hits.k9 <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/ptimefits_plots_top_hits/fits_with_annotations.", jctype, ".k9me3.2022-02-06.neg_slope.txt")
# inf.hits.k27 <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/ptimefits_plots_top_hits/fits_with_annotations.", jctype, ".k27me3.2022-02-05.neg_slope.txt")

for (jctype in jctypes){
  
  
  
}



# Plot bins ---------------------------------------------------------------








