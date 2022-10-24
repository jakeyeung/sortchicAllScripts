# Jake Yeung
# Date of Creation: 2022-07-27
# File: ~/projects/scchic/scripts/revision_scripts/revisions2/9-filter_progenitors_count_tables.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksold) <- jmarks


ctypes.ordered <- c("Granulocytes", "Bcells", "Eryths", "NKs", "Monocytes", "DCs", "Basophils", "MEP", "CMP", "pDCs", "MPPs", "LT", "ST", "HSCs")


# Metas -------------------------------------------------------------------


inf.colors.fixed <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/ctypes_on_umap_batch_corrected_colors_fixed/dat_colors_DC_monocyte_fixed.2022-05-17.txt"
dat.colors.fixed <- fread(inf.colors.fixed)


dat.meta.lst <- lapply(jmarks, function(jmark){
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umaps_pcas_with_batch_corrections/umap_metadata_primetime.", jmark, ".2022-04-21.txt")
  dat.meta <- fread(inf.meta) %>%
    left_join(., dat.colors.fixed) %>%
    rowwise() %>%
    mutate(colcode = colcodenew)
  # replace colcode with colcodenew
})

dat.meta.colors <- subset(dat.meta.lst$k4me1, select = c(ctype.from.LL, colcode))
ctype2col <- hash::hash(dat.meta.colors$ctype.from.LL, dat.meta.colors$colcode)

dat.meta.merge <- lapply(jmarks, function(jmark){
  jdat <- dat.meta.lst[[jmark]]
  subset(jdat, select = c(cell, ctype.from.LL, colcode)) %>%
    rowwise() %>%
    mutate(mark = jmark)
}) %>%
  bind_rows()



# Load dynamic bins -------------------------------------------------------

indir.dynbins <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/first_submission_data/dynamic_bins"

dynbins.lst <- lapply(jmarksold, function(jmarkold){
  print(jmarkold)
  inf.dynbin <- file.path(indir.dynbins, paste0("dynamic_bins_50kb.", jmarkold, ".txt"))
  fread(inf.dynbin)$V4
})

print(lapply(dynbins.lst, length))
old2new <- hash::hash(jmarksold, jmarks)


# Keep progenitors only ---------------------------------------------------

prog.ctypes <- c("CMP", "MPPs", "LT", "ST", "HSCs")
cells.keep.lst <- lapply(jmarks, function(jmark){
  cells.keep <- subset(dat.meta.lst[[jmark]], ctype.from.LL %in% prog.ctypes & batch == "New")$cell
})

lapply(cells.keep.lst, length)

# Load count mats ---------------------------------------------------------

inmain.mat <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/count_mats_bins/counts_tables_50000"
count.mat.dynbins.progs.lst <- lapply(jmarksold, function(jmarkold){
  print(jmarkold)
  indir.mat <- file.path(inmain.mat, paste0("BM_", jmarkold))
  inf.mat <- file.path(indir.mat, paste0("BM_allmerged_", jmarkold, ".countTable.binsize_50000.csv"))
  mat.allbins <- scchicFuncs::ReadMatSlideWinFormat(inf.mat, as.sparse = TRUE, sort.rnames = TRUE, add.chromo = TRUE)
  print(dim(mat.allbins))
  jmarknew <- old2new[[jmarkold]]
  mat.dynbins <- mat.allbins[dynbins.lst[[jmarknew]], ]
  print(dim(mat.dynbins))
  mat.dynbins.progs <- mat.dynbins[, cells.keep.lst[[jmarknew]]]
  print(dim(mat.dynbins.progs))
  return(mat.dynbins.progs)
})



# Save outputs ------------------------------------------------------------

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/count_mats_progenitors_filt"

for (jmark in jmarks){
  print(jmark)
  fout <- file.path(outdir, paste0("countmat_dynbins_new_progenitors_only.", jmark, ".", Sys.Date(), ".rds"))
  saveRDS(count.mat.dynbins.progs.lst[[jmark]], file = fout)
}

print(lapply(count.mat.dynbins.progs.lst, dim))


# Do TSS ------------------------------------------------------------------

matmain.tss <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/allmerged/count_mats_TSS/counts_tables_10000"
count.mat.tss.progs.lst <- lapply(jmarksold, function(jmarkold){
  print(jmarkold)
  matdir <- file.path(matmain.tss, paste0("BM_", jmarkold))
  inf <- file.path(matdir, paste0("BM_allmerged_", jmarkold, ".countTable.binsize_10000.csv"))
  print(inf)
  mat.tss <- ReadMatTSSFormat(inf = inf, as.sparse = TRUE, add.coord = TRUE, sort.rnames = TRUE)
  print(dim(mat.tss))
  jmarknew <- old2new[[jmarkold]]
  mat.tss.progs <- mat.tss[, cells.keep.lst[[jmarknew]]]
  print(dim(mat.tss.progs))
  return(mat.tss.progs)
})

print(lapply(count.mat.tss.progs.lst, dim))

for (jmark in jmarks){
  print(jmark)
  fout <- file.path(outdir, paste0("countmat_TSS_new_progenitors_only.", jmark, ".", Sys.Date(), ".rds"))
  saveRDS(count.mat.tss.progs.lst[[jmark]], file = fout)
}



