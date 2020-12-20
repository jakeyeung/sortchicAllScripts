# Jake Yeung
# Date of Creation: 2020-11-23
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/13-write_glmpca_factors_for_scot.R
# Write latent variables for calculating distances later

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)
library(scchicFuncs)

library(hash)
library(igraph)
library(umap)

library(glmpca)


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")

# Load LDA output ---------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

# niter <- "1000"
# binskeep <- 0
niter <- "500"
binskeep <- 1000
jsuffix <- paste0("bincutoff_0.binskeep_", binskeep, ".byplate.szname_none.niter_", niter, ".reorder_rownames.dupfilt")
# outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2"
# outname <- paste0("bonemarrow_celltypes.", Sys.Date(), ".niter_", niter, ".pdf")
# outname <- paste0("BM_celltypes.", jsuffix, ".", Sys.Date(), ".pdf")
# outf <- file.path(outdir, outname)

infs.lst <- lapply(jmarks, function(jmark){
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/glmpca.", jmark, ".", jsuffix, ".RData"))
  assertthat::assert_that(file.exists(inf))
  return(inf)
})


glmpca.factors.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- infs.lst[[jmark]]
  # inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/glmpca_outputs/glmpca.", jmark, ".", jsuffix, ".RData"))
  load(inf, v=T)
  return(glm.out$factors)
})


# Write to output ---------------------------------------------------------

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_outputs_factors_for_integration"
for (jmark in jmarks){
  print(jmark)
  fname <- paste0("glmpca_factors_", jmark, ".niter_", niter, ".binskeep", binskeep, ".txt")
  fwrite(glmpca.factors.lst[[jmark]], file = file.path(outdir, fname), quote = FALSE, sep = "\t")
}

# write rownames
for (jmark in jmarks){
  print(jmark)
  fname <- paste0("glmpca_factors_", jmark, ".niter_", niter, ".binskeep", binskeep, ".rownames.txt")
  outf <- file.path(outdir, fname)
  assertthat::assert_that(!file.exists(outf))
  fwrite(as.data.frame(rownames(glmpca.factors.lst[[jmark]])), file = outf, col.names = FALSE)
}


# Write celltypes ---------------------------------------------------------

# # update UMAPs
# indir.umaps <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2"
# dat.umaps <- lapply(jmarks, function(jmark){
#   fname.tmp <- paste0("BM_celltypes.bincutoff_0.binskeep_1000.byplate.szname_none.niter_500.reorder_rownames.dupfilt.2020-11-23.", jmark, ".txt")
#   inf.umaps <- file.path(indir.umaps, fname.tmp)
#   dat.tmp <- fread(inf.umaps)
#   dat.tmp$mark <- jmark
#   dat.tmp <- subset(dat.tmp, select = c(cell, umap1, umap2))
#   return(dat.tmp)
# })


