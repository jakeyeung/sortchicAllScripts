# Jake Yeung
# Date of Creation: 2022-02-08
# File: ~/projects/scchic/scripts/scripts_analysis/revisions_from_istbea/12-combine_mats_from_TSS_k27me3.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

jmarks <- c("k27me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K27me3"); names(jmarksold) <- jmarks
outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/count_mat_TSS_rep2rep3only"
inmain <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables_fromTSS/counts_tables_10000"

for (jmark in jmarks){
  
  jmarkold <- jmarksold[[jmark]]
  
  
  indir <- file.path(inmain, paste0("BM_", jmark))
  
  # Load metas --------------------------------------------------------------
  
  # indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/multinom_celltyping/rename"
  # dat.meta <- fread(file.path(indir.meta, paste0("metadata_", jmark, ".txt")))
  # indir.meta <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables/filtered_count_tables_for_LDA_no_frip_filter_on_eryths/BM_k27me3"
  # dat.meta <- fread(file.path(indir.meta, paste0("qc_metadata_new_only.", jmark, ".2022-04-12.txt"))) %>%
  #   filter(is.good)
  # 
  # # load old mat
  # cells.keep.new <- dat.meta$cell
  mat.with.good.cells <- readRDS("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/raw_data/count_tables/filtered_count_tables_for_LDA_no_frip_filter_on_eryths/BM_k27me3/count_mat_merged_with_old_dynbins.k27me3.2022-04-12.rds")
  cells.keep <- colnames(mat.with.good.cells)
  
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
  
  # inf.lda <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.from_TSS.dist_10000.RemoveDupRows/lda_outputs.count_mat_from_TSS.", jmarkold, ".dist_10000.K-30.binarize.FALSE/ldaOut.count_mat_from_TSS.", jmarkold, ".dist_10000.K-30.Robj")
  inf.lda <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_TSS/lda_outputs.PZ-BM-rep3-", jmarkold, "-rep2rep3reseq.TSS.varfilt.K-30.binarize.FALSE/ldaOut.PZ-BM-rep3-", jmarkold, "-rep2rep3reseq.TSS.varfilt.K-30.Robj")
  load(inf.lda, v=T)
  count.mat.old <- count.mat
  
  # Concatenate rows -------------------------------------------------------------
  
  # annotate rnames
  rnames.annotate <- sapply(rownames(count.mat.new), function(x) paste0("chr", strsplit(x, ";")[[1]][[1]]))
  rnames.hash <- hash::hash(rnames.annotate, names(rnames.annotate))
  
  # rename countmat
  rownames(count.mat.old) <- sapply(rownames(count.mat.old), function(x) AssignHash(x = x, jhash = rnames.hash, null.fill = x))
  
  rnames.common <- intersect(rownames(count.mat.new), rownames(count.mat.old))
  
  print(length(rnames.common))
  
  count.mat.merged <- cbind(count.mat.old[rnames.common, ], count.mat.new[rnames.common, ])
  # cells.keep.i <- colnames(count.mat.merged) %in% cells.keep
  count.mat.merged.filt <- count.mat.merged[, cells.keep]
  
  print(dim(count.mat.merged.filt))
  print(length(cells.keep))
  
  # save outputs
  outrds <- file.path(outdir, paste0("count_mat_TSS_rep2rep3only.", jmark, ".", Sys.Date(), ".rds"))
  saveRDS(object = count.mat.merged.filt, file = outrds)
  
  # Filter out cells  -------------------------------------------------------
  
}




