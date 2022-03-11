# Jake Yeung
# Date of Creation: 2022-02-07
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/12-combine_mats_from_TSS.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksold) <- jmarks
outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mat_TSS"
inmain <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables_fromTSS/counts_tables_10000"

for (jmark in jmarks){
  
  jmarkold <- jmarksold[[jmark]]
  
  
  indir <- file.path(inmain, paste0("BM_", jmark))
  
  # Load metas --------------------------------------------------------------
  
  indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping/rename"
  dat.meta <- fread(file.path(indir.meta, paste0("metadata_", jmark, ".txt")))
  
  cells.keep <- dat.meta$cell
  
  # Load mats from new ------------------------------------------------------
  
  infs.mat <- list.files(indir, pattern = ".csv", full.names = TRUE)
  
  library(scchicFuncs)
  mats.lst <- lapply(infs.mat, function(jinf){
    m <- ReadMatTSSFormat(jinf, add.coord = TRUE)
    return(m)
  })
  
  rnames.keep <- gtools::mixedsort(unique(unlist(sapply(mats.lst, function(x) rownames(x)))))
  
  count.mat.new <- cbind.fill.lst(mats.lst, all.rnames = rnames.keep, fill = 0)
  
  # Load mats from old ------------------------------------------------------
  
  inf.lda <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.from_TSS.dist_10000.RemoveDupRows/lda_outputs.count_mat_from_TSS.", jmarkold, ".dist_10000.K-30.binarize.FALSE/ldaOut.count_mat_from_TSS.", jmarkold, ".dist_10000.K-30.Robj")
  load(inf.lda, v=T)
  count.mat.old <- count.mat
  
  
  # Concatenate rows -------------------------------------------------------------
  
  rnames.common <- intersect(rownames(count.mat.new), rownames(count.mat.old))
  
  print(length(rnames.common))
  
  count.mat.merged <- cbind(count.mat.old[rnames.common, ], count.mat.new[rnames.common, ])
  cells.keep.i <- colnames(count.mat.merged) %in% cells.keep
  count.mat.merged.filt <- count.mat.merged[, cells.keep.i]
  
  print(dim(count.mat.merged.filt))
  print(length(cells.keep))
  
  # save outputs
  outrds <- file.path(outdir, paste0("count_mat_TSS_combined.", jmark, ".", Sys.Date(), ".rds"))
  saveRDS(object = count.mat.merged.filt, file = outrds)
  
  # Filter out cells  -------------------------------------------------------
  
  
  
  
}




