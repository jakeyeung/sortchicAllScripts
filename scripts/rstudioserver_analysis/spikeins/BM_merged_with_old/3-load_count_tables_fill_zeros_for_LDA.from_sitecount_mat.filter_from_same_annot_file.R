# Jake Yeung
# Date of Creation: 2020-11-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/3-load_count_tables_fill_zeros_for_LDA.from_sitecount_mat.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

cbind.fill<-function(..., fill = 0){
  nm <- list(...) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(data = fill, n-nrow(x), ncol(x))))) 
}

cbind.fill.lst <- function(mats.lst, all.rnames, fill = 0){
  mats.lst.filled <- lapply(mats.lst, function(mat.tmp){
    missing.rnames <- all.rnames[!all.rnames %in% rownames(mat.tmp)]
    mat.tmp.to.fill <- matrix(data = fill, nrow = length(missing.rnames), ncol = ncol(mat.tmp), dimnames = list(missing.rnames, colnames(mat.tmp)))
    mat.tmp.bind <- rbind(mat.tmp, mat.tmp.to.fill)
    mat.tmp.bind <- mat.tmp.bind[all.rnames, ]
    return(mat.tmp.bind)
  })
  return(do.call(cbind, mats.lst.filled))
}

hubprefix <- "/home/jyeung/hub_oudenaarden"
outdir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BMround2.from_peaks.sitecount_mat")
dir.create(outdir)

outdir2 <- file.path(outdir, "filtNAcells_allbins.from_same_annot_file")
dir.create(outdir2)

# indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second/count_tables_from_peaks")
indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/count_tables_from_peaks.from_sitecount_mat")
assertthat::assert_that(dir.exists(indir))

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

mats.bymark.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  infs.bymark <- list.files(path = indir, pattern = paste0("BM_round1_round2_merged_", jmark, ".*.txt"), full.names = TRUE)
  names(infs.bymark) <- sapply(infs.bymark, basename)
  
  mats.bymark <- lapply(infs.bymark, function(inf){
    ReadMatTSSFormat(inf, add.coord = TRUE, sort.rnames = TRUE)
  })
  
  jall.rnames <- sort(unique(unlist(lapply(mats.bymark, rownames))))
  mats.bymark.fill <- cbind.fill.lst(mats.bymark, all.rnames = jall.rnames, fill = 0)
  return(mats.bymark.fill)
})


# Write to output ---------------------------------------------------------

for (jmark in jmarks){
  print(jmark)
  outf <- file.path(outdir, paste0("count_mat_from_sitecount_mat.", jmark, ".rds"))
  saveRDS(mats.bymark.lst[[jmark]], file = outf)
}


# Filter NA cells and rewrite ---------------------------------------------


# jtest <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins/cell_cluster_table_with_spikeins.H3K4me1.2020-11-18.dupfilt.txt"
# jdat <- fread(jtest)
# subset(jdat, cluster == "Basophils")


# cells.dir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BMAllMerged2.from_peaks/filtNAcells_allbins")
cells.dir <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins")

cells.keep.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- file.path(cells.dir, paste0("cell_cluster_table_with_spikeins.", jmark, ".2020-11-18.dupfilt.txt"))
  # fname <- paste0("count_mat_from_hiddendomains.", jmark, ".filtNAcells_allbins.rds")
  dat.check <- fread(fname)
  # dat.check <- readRDS(file.path(cells.dir, fname))
  cells.keep <- dat.check$cell
  # cells.keep <- colnames(dat.check)
})


for (jmark in jmarks){
  print(jmark)
  cells.keep <- cells.keep.lst[[jmark]]
  outf2 <- file.path(outdir2, paste0("count_mat_from_sitecount_mat.", jmark, ".filtNAcells_allbins.from_same_annot_file.rds"))
  print(dim(mats.bymark.lst[[jmark]]))
  jmat.tmp <- mats.bymark.lst[[jmark]][, cells.keep]
  print(dim(jmat.tmp))
  saveRDS(jmat.tmp, file = outf2)
}



