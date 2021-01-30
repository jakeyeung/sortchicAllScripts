# Jake Yeung
# Date of Creation: 2021-01-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_10kb_analysis/1-load_mat_save_rds.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

indir.cuts.other <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/count_tables_from_bins/binsize_10000")
indir.cuts.k27me3 <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/bams_split_by_cluster/bams_merged_by_cluster/counts_tables_10000")


# load raw cuts

count.mat.lst <- lapply(jmarks, function(jmark){
  indir.cuts <- ifelse(jmark == "H3K27me3", indir.cuts.k27me3, indir.cuts.other)
  assertthat::assert_that(dir.exists(indir.cuts))
  infs.cuts <- list.files(path = indir.cuts, pattern = paste0(jmark, ".*.txt"), all.files = TRUE, full.names = TRUE)
  
  dat.csv.lst <- lapply(infs.cuts, function(inf.tmp){
    print(inf.tmp)
    # mat.tmp <- ReadMatTSSFormat(inf.tmp, as.sparse = TRUE, add.coord = TRUE, sort.rnames = FALSE)
    mat.tmp <- ReadMatSlideWinFormat(inf.tmp, as.sparse = TRUE, sort.rnames = FALSE, add.chromo = TRUE)
    rownames(mat.tmp) <- paste("chr", rownames(mat.tmp), sep = "")
    return(mat.tmp)
  })
  
  rnames.all <- sort(unique(unlist(lapply(dat.csv.lst, function(jmat) rownames(jmat)))))
  count.mat <- cbind.fill.lst(dat.csv.lst, all.rnames = rnames.all)
  
})


print(lapply(count.mat.lst, function(x) x[1:5, 1:5]))

# Filter cells  -----------------------------------------------------------

indir.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins.2020-12-22.umap_spread.H3K27me3_cleaned"
dat.metas <- lapply(jmarks, function(jmark){
  print(jmark)
  fname.tmp <- paste0("cell_cluster_table_with_spikeins.", jmark, ".2020-12-27.umap_spread.final.order_by_cuts_to_spikeins.txt")
  inf <- file.path(indir.meta, fname.tmp)
  jdat <- fread(inf)
  if (jmark == "H3K9me3"){
    jdat$mark <- jmark
  }
  return(jdat)
})


# Write outputs -----------------------------------------------------------

count.mat.filt.lst <- lapply(jmarks, function(jmark){
  cells.keep <- dat.metas[[jmark]]$cell
  cells.keep2 <- colnames(count.mat.lst[[jmark]]) %in% cells.keep
  count.mat.lst[[jmark]][, cells.keep2]
})

lapply(count.mat.lst, dim)
lapply(count.mat.filt.lst, dim)

# save RDS
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_tables.BMfinal.from_bins.10kb"
for (jmark in jmarks){
  print(jmark)
  outf <- file.path(outdir, paste0("count_table_bins.dist_10000.", jmark, ".rds"))
  saveRDS(count.mat.lst[[jmark]], file = outf)
}



