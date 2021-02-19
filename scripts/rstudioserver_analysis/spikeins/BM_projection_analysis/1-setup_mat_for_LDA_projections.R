# Jake Yeung
# Date of Creation: 2020-12-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_projection_analysis/1-setup_mat_for_LDA_projections.R
# Load raw matrices, filter cells, run LDA for rep2rep3 only, but save the rep1old for separate matrix, ready to project.


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

# H3K9me3 is not same_annot_file
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"


# Load meta  --------------------------------------------------------------

dat.metas <- lapply(jmarks, function(jmark){
  if (jmark != "H3K27me3"){
    inf.meta <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned/cell_cluster_table_with_spikeins.", jmark, ".2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt"))
  } else {
    inf.meta <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq.with_old/mat_H3K27me3_rep1rep2rep3reseq.metadata.txt"))
  }
  fread(inf.meta)
})

# Load LDAs ---------------------------------------------------------------


inf.lst <- lapply(jmarks, function(jmark){
  if (jmark %in% c("H3K4me1", "H3K4me3")){
    inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins.from_sitecount_mat.from_same_annot_file/lda_outputs.count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.K-30.binarize.FALSE/ldaOut.count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.K-30.Robj"))
  } else if (jmark == "H3K27me3"){
    # inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_varfilt_peaks/lda_outputs.PZ-BM-rep2rep3reseq-", jmark, ".peaks.varfilt.K-30.binarize.FALSE/ldaOut.PZ-BM-rep2rep3reseq-", jmark, ".peaks.varfilt.K-30.Robj"))
    inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep1rep2rep3reseq_varfilt_2020-12-15/lda_outputs.mat_", jmark, "_rep1rep2rep3reseq.peaks.K-30.binarize.FALSE/ldaOut.mat_", jmark, "_rep1rep2rep3reseq.peaks.K-30.Robj"))
  } else if (jmark == "H3K9me3"){
    inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.frompeaks.filtNAcells_allbins/lda_outputs.count_mat_from_hiddendomains.", jmark, ".filtNAcells_allbins.K-30.binarize.FALSE/ldaOut.count_mat_from_hiddendomains.", jmark, ".filtNAcells_allbins.K-30.Robj"))
  } else {
    print(paste("jmark not found:", jmark))
  }
  print(inf)
  assertthat::assert_that(file.exists(inf))
  return(inf)
})


lda.lst <- lapply(inf.lst, function(inf){
  load(inf, v=T)
  return(list(out.lda = out.lda, count.mat = count.mat))
})


# Filter by rep1old or not ------------------------------------------------

cells.old.lst <- lapply(jmarks, function(jmark){
  subset(dat.metas[[jmark]], jrep == "rep1old")$cell
})

cells.new.lst <- lapply(jmarks, function(jmark){
  subset(dat.metas[[jmark]], jrep != "rep1old")$cell
})

# Write raw counts by matrix ----------------------------------------------

# write for old cells, skip K27me3

outdir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BMround2.from_peaks.sitecount_mat.split_old_and_new")

for (jmark in jmarks){
  print(jmark)
  if (length(cells.old.lst[[jmark]]) == 0){
    print(paste("No cells. Skipping", jmark))
    next
  }
  count.mat.old <- lda.lst[[jmark]]$count.mat
  cnames.old <- colnames(count.mat.old) %in% cells.old.lst[[jmark]]
  count.mat.old <- count.mat.old[, cnames.old]
  print(dim(count.mat.old))
  fnameold <- paste0("count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.old_cells_only.", Sys.Date(), ".rds")
  outoldtxt <- file.path(outdir, fnameold)
  saveRDS(count.mat.old, file = outoldtxt)
}

for (jmark in jmarks){
  print(jmark)
  if (length(cells.new.lst[[jmark]]) == 0){
    print(paste("No cells. Skipping", jmark))
    next
  }
  count.mat.new <- lda.lst[[jmark]]$count.mat
  cnames.new <- colnames(count.mat.new) %in% cells.new.lst[[jmark]]
  count.mat.new <- count.mat.new[, cnames.new]
  print(dim(count.mat.new))
  fnamenew <- paste0("count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.new_cells_only.", Sys.Date(), ".rds")
  outnewtxt <- file.path(outdir, fnamenew)
  saveRDS(count.mat.new, file = outnewtxt)
}


# Write metadata ----------------------------------------------------------

for (jmark in jmarks){
  fname <- paste0("count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.metadata.", Sys.Date(), ".txt")
  outmeta <- file.path(outdir, fname)
  print(jmark)
  fwrite(dat.metas[[jmark]], file = outmeta, sep = "\t")
}





