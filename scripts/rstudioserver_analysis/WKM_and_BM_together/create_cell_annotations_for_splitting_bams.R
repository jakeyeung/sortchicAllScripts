# Jake Yeung
# Date of Creation: 2020-06-17
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/create_cell_annotations_for_splitting_bams.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/cell_to_cluster_tables_merged"

# Load BM -----------------------------------------------------------------


fewer.k27me3 <- TRUE
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.setup_for_poisson_regression"
dir.create(indir, showWarnings = TRUE)
jprefix <- file.path(indir, paste0("integrated_analysis_3_marks.stringentDE.WithHighLowExprs.singlecells.fewerk27me3_", fewer.k27me3, ".forPoissonRegression.CountR1only.", "2020-06-05"))

infrdata <- paste0(jprefix, ".smaller.RData")

load(infrdata, v=T)

dat.annots.filt.BM <- dat.annots.filt

dat.annots.BM.output <- lapply(dat.annots.filt.BM, function(jdat){
  subset(jdat, select = c(cell, cluster.new)) %>%
    dplyr::rename(cluster = cluster.new)
})


# Load ZF  ----------------------------------------------------------------



jdate <- "2020-06-09"
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/setup_matrix_for_poisson_regression.likeBM.redo_count_tables"
assertthat::assert_that(dir.exists(indir))
# dir.create(indir, showWarnings = TRUE)
jprefix <- file.path(indir, paste0("integrated_analysis.", jdate, ".UseTSSfromH3K4me3.likeBM"))
infrdata <- paste0(jprefix, ".RData")
load(infrdata, v=T)

dat.annots.filt.WKM <- dat.annot.lst.WKM

dat.annots.WKM.output <- lapply(dat.annots.filt.WKM, function(jdat){
  subset(jdat, select = c(cell, cluster.new)) %>%
    dplyr::rename(cluster = cluster.new)
})




# Get K9me3 and write -----------------------------------------------------

jmarks.rep <- c("H3K9me3")
names(jmarks.rep) <- jmarks.rep

indir.annots <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping")


dat.annots.bm.k9 <- lapply(jmarks.rep, function(jmark){
  inf.annots <- file.path(indir.annots, paste0("GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData"))
  load(inf.annots, v=T)
  dat.output <- dat.umap.glm.fillNAs %>%
    dplyr::select(cell, cluster)
  return(dat.output)
})

dat.annots.wkm.k9 <- lapply(jmarks.rep, function(jmark){
  print(jmark)
  # filter by previously defined celltypes? 
  inf.annot.louv <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/LDA_downstream/LDA_downstream_ZF.2020-04-23.imputevarfilt.lessstringent/ZF_LDA_output.", jmark, ".keepn_150.final.ClusterTables.txt")
  assertthat::assert_that(file.exists(inf.annot.louv))
  
  inf.annot.glmpca <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/celltyping/from_louvain/cell_to_cluster_table.", jmark, ".2020-04-13.txt")
  assertthat::assert_that(file.exists(inf.annot.glmpca))
  
  annot.louv <- fread(inf.annot.louv)
  annot.louv$clusterplate <- paste(annot.louv$cluster, annot.louv$plate, "_")
  
  annot.glmpca <- fread(inf.annot.glmpca)
  annot.glmpca.filt <- subset(annot.glmpca, cell %in% annot.louv$cell)
  # rename clusters
  annot.glmpca.filt$cond <- sapply(annot.glmpca.filt$cell, ClipLast, jsep = "_")
  annot.output <- subset(annot.glmpca.filt, select = c(cell, cluster))
  return(annot.output)
})

dat.annots.BM.output[[jmarks.rep]] <- dat.annots.bm.k9[[jmarks.rep]]
dat.annots.WKM.output[[jmarks.rep]] <- dat.annots.wkm.k9[[jmarks.rep]]

jmarks.all <- c(jmarks, jmarks.rep)

# Write to outputs --------------------------------------------------------

outprefix.bm <- "BM_cell_to_clusters"
outprefix.wkm <- "WKM_cell_to_clusters"

lapply(jmarks.all, function(jmark){
  outf.tmp.bm <- file.path(outdir, paste0(outprefix.bm, ".", jmark, ".txt"))
  outf.tmp.wkm <- file.path(outdir, paste0(outprefix.wkm, ".", jmark, ".txt"))
  fwrite(dat.annots.BM.output[[jmark]], file = outf.tmp.bm, sep = "\t", na = "NA")
  fwrite(dat.annots.WKM.output[[jmark]], file = outf.tmp.wkm, sep = "\t", na = "NA")
})

