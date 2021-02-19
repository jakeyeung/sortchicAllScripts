# Jake Yeung
# Date of Creation: 2021-01-31
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/16-load_raw_counts_filter_cells_for_LDA_k4_k9_all_50kb_bins.remove_HSPCs.R
# k

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)


GetCountMat <- function(infs.k4){
  # (infs.k4 <- list.files(indir, pattern = "BM_round1_round2_merged_H3K4me1_.*..bam.count_table_k4_k9_dynamic_regions.txt", full.names = TRUE))
  counts.mat.k4.lst <- lapply(infs.k4, function(inf){
    count.mat.tmp <- ReadMatTSSFormat(inf, as.sparse = TRUE, add.coord = TRUE, sort.rnames = FALSE)
    return(count.mat.tmp)
  })
  rnames.all.k4 <- unique(unlist(lapply(counts.mat.k4.lst, function(x) rownames(x))))
  count.mat.k4 <- cbind.fill.lst(counts.mat.k4.lst, all.rnames = rnames.all, fill = 0)
}


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"



# Load metas --------------------------------------------------------------

indir.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage")
assertthat::assert_that(dir.exists(indir.meta))
dat.metas <- lapply(jmarks, function(jmark){
  print(jmark)
  # fname.tmp <- paste0("cell_cluster_table_with_spikeins.", jmark, ".2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt")
  fname.tmp <- paste0("metadata_batch_corrected.arranged_by_lineage.", jmark, ".2021-01-02.txt")
  inf <- file.path(indir.meta, fname.tmp)
  print(inf)
  jdat <- fread(inf)
  if (jmark == "H3K9me3"){
    jdat$mark <- jmark
  }
  return(jdat)
})



# Get dbl cells -----------------------------------------------------------

jsuffix <- "bam.count_table_k4_k9_50kb_genomewide.txt"
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/BM_H3K4me1-H3K9me3_tech_rep_merged/raw_demultiplexed/count_tables_from_H3K4me1_H3K9me3_50kb_genomewide"
(infs.dbl <- list.files(indir, pattern = paste0("PZ-BM-rep3-H3K9me3-H3K4me1-.*.", jsuffix), full.names = TRUE))

# "PZ-BM-rep3-H3K9me3-H3K4me1-1.sorted.tagged.bam.count_table_k4_k9_dynamic_bins_from_heatmap.txt"
assertthat::assert_that(length(infs.dbl) > 0)

counts.mat.dbl.lst <- lapply(infs.dbl, function(inf){
  count.mat.tmp <- ReadMatTSSFormat(inf, as.sparse = TRUE, add.coord = TRUE, sort.rnames = FALSE)
  return(count.mat.tmp)
})

rnames.all <- unique(unlist(lapply(counts.mat.dbl.lst, function(x) rownames(x))))

count.mat.dbl <- cbind.fill.lst(counts.mat.dbl.lst, all.rnames = rnames.all, fill = 0)

# get good cells 
inf <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BM_dbl_reseq.varfilt/lda_outputs.count_mat_cleaned_reseq.H3K4me1_H3K9me3.varfilt.K-30.binarize.FALSE/ldaOut.count_mat_cleaned_reseq.H3K4me1_H3K9me3.varfilt.K-30.Robj")
load(inf, v=T)
cells.keep.dbl <- colnames(count.mat)


count.mat.dbl.filt <- count.mat.dbl[, cells.keep.dbl]


# Get singles -------------------------------------------------------------


infs.k4 <- list.files(indir, pattern = paste0("BM_round1_round2_merged_H3K4me1_.*.", jsuffix), full.names = TRUE)
assertthat::assert_that(length(infs.k4) > 0)
count.mat.k4 <- GetCountMat(infs.k4)

infs.k9 <- list.files(indir, pattern = paste0("BM_round1_round2_merged_H3K9me3_.*.", jsuffix), full.names = TRUE)
assertthat::assert_that(length(infs.k9) > 0)
count.mat.k9 <- GetCountMat(infs.k9)


# Get count mat lst -------------------------------------------------------

count.mat.lst <- list("H3K4me1-H3K9me3" = count.mat.dbl.filt,
                      "H3K4me1" = count.mat.k4,
                      "H3K9me3" = count.mat.k9)

jnames <- names(count.mat.lst)

# Write outputs -----------------------------------------------------------

# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_bins_from_heatmap"
outdir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/from_genomewide.RemoveHSPCs")
dir.create(outdir)


# Remove HSPCs ------------------------------------------------------------

cells.keep.k4 <- subset(dat.metas$H3K4me1, batch != "StemCell")$cell
cells.keep.k9 <- subset(dat.metas$H3K9me3, batch != "StemCell")$cell

lapply(count.mat.lst, dim)

count.mat.lst$H3K4me1 <- count.mat.lst$H3K4me1[, cells.keep.k4]
count.mat.lst$H3K9me3 <- count.mat.lst$H3K9me3[, cells.keep.k9]

lapply(count.mat.lst, dim)

for (jname in jnames){
  print(jname)
  fname <- paste0("count_name.", jname, ".k4_k9_50kb_genomewide.RemoveHSPCs", Sys.Date(), ".rds")
  outf <- file.path(outdir, fname)
  jtmp <- count.mat.lst[[jname]]
  print(dim(jtmp))
  saveRDS(jtmp, file = outf)
}
  
