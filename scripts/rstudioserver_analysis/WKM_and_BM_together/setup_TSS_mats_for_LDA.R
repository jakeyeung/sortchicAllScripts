# Jake Yeung
# Date of Creation: 2020-06-22
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/setup_TSS_mats_for_LDA.R
# Set it up 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/mats_for_LDA"

# Mouse BM -----------------------------------------------------------------------


fewer.k27me3 <- TRUE
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.setup_for_poisson_regression"
dir.create(indir, showWarnings = TRUE)
jprefix <- file.path(indir, paste0("integrated_analysis_3_marks.stringentDE.WithHighLowExprs.singlecells.fewerk27me3_", fewer.k27me3, ".forPoissonRegression.CountR1only.", "2020-06-05"))

infrdata <- paste0(jprefix, ".smaller.RData")

load(infrdata, v=T)

tss.keep.lst <- lapply(tss.mats.filt.fromref.cellfilt, function(x) rownames(x))

lapply(tss.keep.lst, length)

tss.keep.common <- Reduce(f = intersect, x = tss.keep.lst)
# genes.keep.common <- sapply(tss.keep.common, function(x) strsplit(x, ";")[[1]][[2]], USE.NAMES = FALSE)
tx.keep.common <- sapply(tss.keep.common, function(rname) paste(strsplit(rname, split = ";")[[1]][2:3], collapse = ";"), USE.NAMES = FALSE)


# load full TSS matrix again (some cells were =removed in TSS MATS)
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/debug/ZellerRawData_B6_All_MergedByMarks_final.count_tables_TSS.all_tx.noR2.Buys"

mats.50kb <- lapply(jmarks, function(jmark){
  fname <- paste0(jmark, ".countTableTSS.mapq_40.TSS_50000.blfiltered.csv")
  inf <- file.path(indir, fname)
  mat <- ReadMatTSSFormat(inf = inf, as.sparse = TRUE, add.coord = FALSE, sort.rnames = TRUE)
  return(mat)
})

mats.10kb <- lapply(jmarks, function(jmark){
  fname <- paste0(jmark, ".countTableTSS.mapq_40.TSS_10000.blfiltered.csv")
  inf <- file.path(indir, fname)
  mat <- ReadMatTSSFormat(inf = inf, as.sparse = TRUE, add.coord = FALSE, sort.rnames = TRUE)
  return(mat)
})

mats.10kb.filt1 <- lapply(mats.10kb, function(jmat){
  rnames.all <- rownames(jmat)
  genes.all <- sapply(rnames.all, function(x) strsplit(x, ";")[[1]][[2]])
  tx.all <- sapply(rnames.all, function(rname) paste(strsplit(rname, split = ";")[[1]][2:3], collapse = ";"), USE.NAMES = FALSE)
  rnames.keep <- tx.all %in% tx.keep.common
  jmat[rnames.keep, ]
})

# get common rows
rnames.10kb.keep <- Reduce(intersect, lapply(mats.10kb.filt1, rownames))

mats.50kb.filt1 <- lapply(mats.50kb, function(jmat){
  rnames.all <- rownames(jmat)
  tx.all <- sapply(rnames.all, function(rname) paste(strsplit(rname, split = ";")[[1]][2:3], collapse = ";"), USE.NAMES = FALSE)
  rnames.keep <- tx.all %in% tx.keep.common
  jmat[rnames.keep, ]
})
rnames.50kb.keep <- Reduce(intersect, lapply(mats.50kb.filt1, rownames))

mats.10kb.filt2 <- lapply(mats.10kb.filt1, function(jmat){
  jmat[rnames.10kb.keep, ]
})

mats.50kb.filt2 <- lapply(mats.50kb.filt1, function(jmat){
  jmat[rnames.50kb.keep, ]
})


lapply(mats.10kb.filt2, dim)
lapply(mats.50kb.filt2, dim)

jdistnames <- c("kb10", "kb50")
jdists <- c(10000L, 50000L)
names(jdists) <- jdistnames
mats.lst.bm <- list(kb10 = mats.10kb.filt2, kb50 = mats.50kb.filt2)
names(mats.lst.bm) <- jdistnames


# Write matrix to output --------------------------------------------------


jspecies <- "MouseBM"

for (jdistname in jdistnames){
  print(jdistname)
  jmat.tmp.lst <- mats.lst.bm[[jdistname]]
  for (jmark in jmarks){
    jmat.tmp <- jmat.tmp.lst[[jmark]]
    outf <- file.path(outdir, paste0(jspecies, ".", jmark, ".TSSdist_", jdists[[jdistname]], ".", Sys.Date(), ".rds"))
    print(outf)
    saveRDS(jmat.tmp, file = outf)
  }
}


# Zebrafish ---------------------------------------------------------------


jspecies <- "ZebrafishWKM"

jdate <- "2020-06-09"
indir.wkm <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/setup_matrix_for_poisson_regression.likeBM.redo_count_tables"
assertthat::assert_that(dir.exists(indir.wkm))
# dir.create(indir, showWarnings = TRUE)
jprefix <- file.path(indir.wkm, paste0("integrated_analysis.", jdate, ".UseTSSfromH3K4me3.likeBM"))
infrdata <- paste0(jprefix, ".RData")
load(infrdata, v=T)


# take transcripts: can be different window siezs
refmark <- "H3K4me3"
rnames <- rownames(tss.mats.sc.filt.zf[[refmark]])
tx.keep <- sapply(rnames, function(rname) paste(strsplit(rname, split = ";")[[1]][2:3], collapse = ";"), USE.NAMES = FALSE)


# Load the new tables ------------------------------------------------

jdists <- c(50000L, 10000L)
names(jdists) <- jdistnames

# must be TSS tables

mats.lst.wkm <- lapply(jdists, function(jdist){
  indir.wkm.newtss <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables_all/count_tables.TSS.CodingOnly.imputevarfilt.lessstringent.mapq_40.winsize_", jdist, ".r1only")
  tss.mats.filt.fromref.cellfilt <- lapply(jmarks, function(jmark){
    print(jmark)
    fname <- paste0(jmark, ".imputevarfilt.lessstringent.mapq_40.remerged.countTable.TSS.csv")
    inf <- file.path(indir.wkm.newtss, fname)
    mat <- ReadMatTSSFormat(inf = inf, as.sparse = TRUE, sort.rnames = TRUE)
    # filter rownames to match tx.keep
    tx.all <- sapply(rownames(mat), function(rname) paste(strsplit(rname, ";")[[1]][2:3], collapse = ";"))
    tx.filt <- tx.all %in% tx.keep
    # all celltypes
    return(mat[tx.filt, ])
  })
})

# filter out common rownames
mats.lst.wkm.filt <- lapply(mats.lst.wkm, function(jmats.bymark){
  rows.common <- Reduce(intersect, lapply(jmats.bymark, rownames))
  jmats.bymark.filt <- lapply(jmats.bymark, function(jmat){
    jmat[rows.common, ]
  })
})

lapply(mats.lst.wkm.filt, function(jmats){
  lapply(jmats, dim)
})



# Write to output ---------------------------------------------------------

for (jdistname in jdistnames){
  jdist <- jdist[[jdistname]]
  print(jdistname)
  jmat.tmp.lst <- mats.lst.wkm.filt[[jdistname]]
  for (jmark in jmarks){
    jmat.tmp <- jmat.tmp.lst[[jmark]]
    outf <- file.path(outdir, paste0(jspecies, ".", jmark, ".TSSdist_", jdist, ".", Sys.Date(), ".rds"))
    print(outf)
    assertthat::assert_that(!file.exists(outf))
    saveRDS(jmat.tmp, file = outf)
  }
}

