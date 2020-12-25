# Jake Yeung
# Date of Creation: 2020-12-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/3-bone_marrow_celltypes.from_glmpca.clean_up_K27me3_batch_effects.make_matrices.downsteram.H3K27me3.R
# description


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

jstart <- Sys.Date()


# Load objs  --------------------------------------------------------------

jmarks <- c("H3K27me3", "H3K4me1", "H3K4me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmark.test <- "H3K27me3"
# jmark.test <- "H3K4me1"

mat.adj.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  # inrdata <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/plate_cluster_adj2.H3K4me1_TSS_topics/matrix_imputed_mark_", jmark, ".RData")
  if (jmark == "H3K27me3"){
    inrdata <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/plate_cluster_adj2.k4me1_and_k4me3_genes.K27me3only/matrix_imputed_mark_H3K27me3.RData")
  } else {
    inrdata <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/plate_cluster_adj2.k4me1_and_k4me3_genes/matrix_imputed_mark_", jmark, ".RData")
  }
  load(inrdata, v=T)
  return(mat.adj)
})


# # inrdata <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/plate_cluster_adj2/matrix_imputed_mark_", jmark.test, ".RData")
# 
# 
# # inf.genes <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/with_imputed/celltype_specific_genes_defined_by_K4me3_TSS.txt")
# # inf.genes <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/")


inf.genes.k4me1 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/geneset_from_H3K4me1/geneset_H3K4me1_TSS_topics2.filt.colorcoded2.txt"
dat.genes.k4me1 <- fread(inf.genes.k4me1)

dat.genes.k4me1$coord <- paste("chr", sapply(dat.genes.k4me1$gene, function(x) strsplit(x, ";")[[1]][[1]]), sep = "")
baso.genes <- subset(dat.genes.k4me1, jset == "Basophils" & rnk < 150)
# inf.genes <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/"
inf.genes <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/with_imputed/celltype_specific_genes_defined_by_K4me3_TSS.txt")
dat.genes <- fread(inf.genes)

dat.meta.lst <- lapply(jmarks, function(jmark){
  if (jmark != "H3K27me3"){
    # inf.meta <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/with_imputed/metadata_umap_celltype_cuts.", jmark, ".txt"))
    inf.meta <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/plate_cluster_adj2.k4me1_and_k4me3_genes.K27me3only/matrix_imputed_mark_H3K27me3.newmetadata.pDCstoBasophils.txt"))
  } else {
    inf.meta <- file.path(hubprefix , "jyeung/data/scChiC/from_rstudioserver/count_tables.BMrep2rep3reseq.with_old/mat_H3K27me3_rep1rep2rep3reseq.metadata.txt")
  }
  dat <- fread(inf.meta)
  return(dat)
})


library(DescTools)



dat.genes.all <- bind_rows(dat.genes, subset(dat.genes.k4me1, select = c(-topic)))

dat.genes.sub <- subset(dat.genes, !(jset == "Basophils" & rnk > 100))
dat.genes.sub <- subset(dat.genes.sub, !(jset == "pDCs" & rnk > 100))

dat.genes.sub <- subset(dat.genes, rnk < 150)

# swap basophils
dat.genes.sub <- subset(dat.genes.sub, jset != "Basophils" & jset != "pDCs")

basos.k4me3 <- subset(dat.genes, jset %in% c("Basophils"))$gene
basos.k4me1 <- subset(dat.genes.k4me1, jset %in% c("Basophils"))$gene

pdcs.k4me3 <- subset(dat.genes, jset %in% c("pDCs"))$gene
pdcs.k4me1 <- subset(dat.genes.k4me1, jset %in% c("pDCs"))$gene

basos.manual <- c("Il4", "Il6", "Cpa3", "Il1r1")
basos.rname <- subset(dat.genes.all, symbol %in% basos.manual)$gene
basos.both <- c(intersect(basos.k4me1, basos.k4me3), basos.rname)
pdcs.both <- intersect(pdcs.k4me1, pdcs.k4me3)

# basos.pdcs <- subset(dat.genes, jset %in5 c("pDCs", "Basophils"))$gene

# dat.to.add <- subset(dat.genes.k4me1, jset %in% c("Basophils", "pDCs") & rnk < 50)
dat.to.add <- subset(dat.genes.all, gene %in% c(basos.both, pdcs.both))

dat.genes.sub.join <- bind_rows(dat.genes.sub, dat.to.add)

dat.genes.sub.join$jset <- factor(dat.genes.sub.join$jset, levels = c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "pDCs", "DCs", "HSPCs"))

dat.genes.sub.join <- dat.genes.sub.join %>%
  arrange(jset, rnk)


# label famou("Erythss genes 
famous.genes.lst <- 
  list("Eryths" = c("Tal1", "Aqp1", "Gata1", "Comt", "Sphk1", "Hbb-bs", "Hbb-bt", "Sox6"),
       "Bcells" = c("Pax5", "Ebf1", "Cd79a", "Cd79b", "Snn", "Blnk", "Cd72", "Blk", "Kzf3", "Cd19"), 
       "NKs" = c("Stat4", "Tcf7", "Tbx21", "Cd3d", "Gimap4"), 
       "Granulocytes" = c("Cxcr2", "Mmp9", "Cd177", "Ncf4", "S100a9", "S100a8", "Ncf1", "Cebpe"),
       "Basophils" = c("Il4", "Il1r1", "Arap3", "Cpa3"), 
       "pDCs" = c("Irf8", "Selpg", "Itgb7", "Ccr9", "Unc93b1"), 
       "DCs" = c("Ccl9", "Apoe", "Nlrp3", "Csf1r"), 
       "HSPCs" = c("Hoxa9", "Hoxa7", "Hoxa3", "Meis1", "Runx2", "Kit", "Hlf", "Hoxa10", "Erg", "Cd34", "Hoxa6", "Gata3", "Hoxb4"))

famous.genes <- unlist(famous.genes.lst)

# label cells
dat.genes.sub.join <- dat.genes.sub.join %>%
  rowwise() %>%
  mutate(genesymbol = ifelse(symbol %in% famous.genes, symbol, NA))
# dat.genes.sub.join$genesymbol <- make.names(dat.genes.sub.join$genesymbol, unique = TRUE)



# pdf(paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/bonemarrow_heatmaps_3_marks_TSS_pretty.", Sys.Date(), ".cexsmall2.pdf"), useDingbats = FALSE)
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/heatmap_pdfs_and_ordered_matrices.H3K27me3_with_old"
dir.create(outdir)
outpdf <- file.path(outdir, paste0("bonemarrow_heatmaps_3_marks_TSS_pretty.", Sys.Date(), ".cexsmall.H3K27me3_with_old.pdf"))

pdf(outpdf, useDingbats = FALSE)

for (jmark.test in jmarks){
  outrdata.tmp <- file.path(outdir, paste0("heatmap_ordered_with_labels.", jmark.test, ".", Sys.Date(), ".RData"))
  print(jmark.test)
  
  dat.meta.sub <- dat.meta.lst[[jmark.test]]
  dat.meta.sub$cluster <- factor(dat.meta.sub$cluster, levels = c("Eryths", "Bcells", "NKs", "Granulocytes", "Basophils", "pDCs", "DCs", "HSPCs"))
  dat.meta.sub <- dat.meta.sub %>%
    arrange(cluster, jrep)
  
 
  cells.keep <- dat.meta.sub$cell
  genes.keep <- dat.genes.sub.join$gene
  
  # mat.adj <- as.data.frame(mat.adj)
  mat.adj <- mat.adj.lst[[jmark.test]]
  if ("rname" %in% colnames(mat.adj)){
    mat.adj <- as.data.frame(mat.adj) 
    rownames(mat.adj) <- mat.adj$rname
    mat.adj$rname <- NULL
  }
  # rownames(mat.adj) <- mat.adj$rname
  # mat.adj$rname <- NULL
  
  # winsorize the genes 
  mat.adj.tmp <- mat.adj[genes.keep, cells.keep]
  # rownames(mat.adj.tmp) <- dat.genes.sub.join$genesymbol
  
  if (jmark.test == "H3K27me3"){
    mat.adj.tmp <- t(apply(mat.adj.tmp, 1, function(jrow) Winsorize(jrow, probs = c(0.01, 0.95))))
    mat.adj.tmp <- apply(mat.adj.tmp, 2, function(jcol) Winsorize(jcol, probs = c(0.01, 0.95)))
  } else {
    mat.adj.tmp <- t(apply(mat.adj.tmp, 1, function(jrow) Winsorize(jrow, probs = c(0.05, 0.99))))
    mat.adj.tmp <- apply(mat.adj.tmp, 2, function(jcol) Winsorize(jcol, probs = c(0.05, 0.99)))
    
  }
  labrows = dat.genes.sub.join$genesymbol
  colsidecolors = dat.meta.sub$colorcode
  rowsidecolors = dat.genes.sub.join$colorcode
  heatmap3::heatmap3(mat.adj.tmp, cexRow = 0.08, labRow = labrows, labCol = "", Rowv = NA, Colv = NA, ColSideColors = colsidecolors, scale = "row", RowSideColors = rowsidecolors, revC = TRUE, main = paste0(jmark.test))
  # write matrix and labels to output
  save(mat.adj.tmp, labrows, colsidecolors, rowsidecolors, file = outrdata.tmp)
}
dev.off()


# Write ordered matrix to output ------------------------------------------


