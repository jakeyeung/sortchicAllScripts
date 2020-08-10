# Jake Yeung
# Date of Creation: 2020-06-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/QuarantineAnalyses/setup_tables_for_HSC_LDA.R
# Look for variability within HSCs

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs"
dir.create(outdir)

# Load cell cluster annots ------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

indir.annots <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping")


dat.annots.all <- lapply(jmarks, function(jmark){
  inf.annots <- file.path(indir.annots, paste0("GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData"))
  load(inf.annots, v=T)
  return(dat.umap.glm.fillNAs)
})


# Filter for HSC (two ways) ----------------------------------

jcond <- "StemCell"
dat.annots.condfilt <- lapply(dat.annots.all, function(jdat){
  subset(jdat, cond == jcond)
})

inf.pseudobulk.annot <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.stringentDE.forAvO/filtered_cells_for_pseudobulk.2020-05-31.rds")
dat.pbulk.annot <- readRDS(inf.pseudobulk.annot)

clstr2cond <- lapply(dat.pbulk.annot, function(jdat){
  unique(subset(jdat, cluster.new == "HSPCs")$cluster)
})
clstr2cond$H3K9me3 <- "HSCs_topic22"

# from Louvain, more stringent
dat.annots.louvfilt <- lapply(jmarks, function(jmark){
  jdat <- dat.annots.all[[jmark]]
  jclst <- clstr2cond[[jmark]]
  return(subset(jdat, cluster %in% jclst))
})


# Cells all ---------------------------------------------------------------

cells.all <- lapply(dat.annots.all, function(jdat){
  jdat$cell
})

cells.condfilt <- lapply(dat.annots.condfilt, function(jdat){
  jdat$cell
})

cells.louvfilt <- lapply(dat.annots.louvfilt, function(jdat){
  jdat$cell
})

# Load raw tables ---------------------------------------------------------

indir.tbls <- file.path(hubprefix, "jyeung/data/scChiC/raw_data/debug/ZellerRawData_B6_All_MergedByMarks_final.count_tables_SlidingWindow.slurm.noR2.Buys")
bsize <- 50000
mats.cfilt <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.mat <- file.path(indir.tbls, paste0(jmark, ".mapq_40.SlidingWindow_", bsize, ".blfiltered.csv"))
  assertthat::assert_that(file.exists(inf.mat))
  jmat <- ReadMatSlideWinFormat(inf.mat, as.sparse = TRUE, sort.rnames = TRUE, add.chromo = FALSE)
  jmat <- jmat[, cells.all[[jmark]]]
  return(jmat)
})

lapply(mats.cfilt, dim)

# Remove bad bins ---------------------------------------------------------

exprs.means <- lapply(mats.cfilt, function(jmat){
  rowMeans(jmat)
})

lapply(jmarks, function(jmark){
  x <- exprs.means[[jmark]]
  plot(density(unlist(x)), main = jmark)
})

maxcounts <- 1.25

bins.keep <- lapply(exprs.means, function(jexprs){
  names(jexprs)[which(jexprs <= maxcounts)]
})

bins.keep.intersect <- Reduce(f = intersect, x = bins.keep)

# Filter bins -------------------------------------------------------------

mats.cb.filt <- lapply(jmarks, function(jmark){
  mats.cfilt[[jmark]][bins.keep.intersect, ]
})

lapply(mats.cb.filt, dim)


# Write tables ready for LDA  ---------------------------------------------

# write all
fname.base <- paste0("count_mat_cbfilt_maxcountsfilt.all.", Sys.Date())
lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0(fname.base, ".", jmark, ".rds")
  outf <- file.path(outdir, fname)
  print(outf)
  saveRDS(object = mats.cb.filt[[jmark]], file = outf)
})

# write cond filt
fname.base <- paste0("count_mat_cbfilt_maxcountsfilt.condfilt.", Sys.Date())
lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0(fname.base, ".", jmark, ".rds")
  outf <- file.path(outdir, fname)
  print(outf)
  jtmp <- mats.cb.filt[[jmark]][, cells.condfilt[[jmark]]]
  print(dim(jtmp))
  saveRDS(object = jtmp, file = outf)
})

# write louvfilt
fname.base <- paste0("count_mat_cbfilt_maxcountsfilt.louvfilt.", Sys.Date())
lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0(fname.base, ".", jmark, ".rds")
  outf <- file.path(outdir, fname)
  jtmp <- mats.cb.filt[[jmark]][, cells.louvfilt[[jmark]]]
  print(dim(jtmp))
  saveRDS(object = jtmp, file = outf)
})





