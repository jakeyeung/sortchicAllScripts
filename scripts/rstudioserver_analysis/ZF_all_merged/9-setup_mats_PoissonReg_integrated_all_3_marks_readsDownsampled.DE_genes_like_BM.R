t# Jake Yeung
# Date of Creation: 2020-05-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/8-finalize_downsampled_TSS_counts_integrated_all_3_marks_readsDownsampled.DE_genes_like_BM.R
# 

rm(list=ls())

jstart <- Sys.time()

set.seed(0)

library(DropletUtils)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(purrr)

library(hash)
library(igraph)
library(umap)

library(JFuncs)
library(scchicFuncs)

library(topicmodels)

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ChIPseeker)
library(GenomicRanges)


library(ggrastr)

# Constants ---------------------------------------------------------------

wkm.rename <- hash(c("eryth1", "eryth2", "HSC1", "HSC2", "monocyte"), c("eryth", "eryth", "HSPCs", "HSPCs", "granu"))

jwinsize <- "10000"


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
jvarcutoffs <- c(0.75, 2, 1, 0.5)
names(jvarcutoffs) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

jchromos <- paste("chr", seq(25), sep = "")

jsize <- 4
cbPalette <- c("grey85", "#32CD32", "#56B4E9", "#CC79A7", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/setup_matrix_for_poisson_regression.likeBM"
dir.create(outdir)

outprefix <- paste0("integrated_analysis.", Sys.Date(), ".UseTSSfromH3K4me3.likeBM.")
outname <- paste0(outprefix, ".pdf")
outname.rdata <- paste0(outprefix, ".RData")

outpdf <- file.path(outdir, outname)
outrdata <- file.path(outdir, outname.rdata)



# Read TSS Signal to figure out which transcript to keep  -----------------

indir.tss <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/count_tables.TSS.CodingOnly.winsize_", jwinsize)
assertthat::assert_that(dir.exists(indir.tss))

tss.out <- lapply(jmarks, function(jmark){
  print(jmark)
  inf.tss <- file.path(indir.tss, paste0("PZ-ChIC-ZF_", jmark, "_2020-04-07.countTable.TSS.csv"))
  mat.tss <- ReadMatTSSFormat(inf.tss)
  # mat.tss <- CollapseRowsByGene(mat.tss, as.long=FALSE, track.kept.gene = TRUE)  # do this for just one mark
  return(list(mat.tss = mat.tss, tss.exprs = rowSums(mat.tss)))
})

tss.exprs.lst.unfilt <- lapply(tss.out, function(x) x$tss.exprs)
tss.mats.singlecell.unfilt <- lapply(tss.out, function(x) x$mat.tss)



# exprs.vec <- tss.exprs.lst$H3K4me1

lapply(jmarks, function(jmark){
  plot(density(tss.exprs.lst.unfilt[[jmark]]), main = jmark)
})


# go with the K4me3 definition...

ref.mark <- "H3K4me3"
jthres <- 50  # maybe not exactly at hump? what about tissuespecific stuff? rare celltypes? complicated from the bulk 
plot(density(tss.exprs.lst.unfilt[[ref.mark]]))
abline(v = jthres)

tss.mat.ref <- CollapseRowsByGene(count.mat = tss.mats.singlecell.unfilt[[ref.mark]], as.long = FALSE, track.kept.gene = TRUE)

tss.keep <- rownames(tss.mat.ref)

# # pick TSS's that are active, the rest are thrown out 
# tss.keep.i <- which(tss.exprs.lst[[ref.mark]] >= jthres)
# tss.keep <- tss.exprs.lst[[ref.mark]][tss.keep.i]

tss.exprs.lst <- lapply(tss.exprs.lst.unfilt, function(exprs.vec){
  jkeep <- names(exprs.vec) %in% tss.keep
  return(exprs.vec[jkeep])
})

print("Dimensions of TSS raw keeping all TSS")
lapply(tss.mats.singlecell.unfilt, length)
tss.mats.singlecell <- lapply(tss.mats.singlecell.unfilt, function(tss.mat){
  jkeep <- rownames(tss.mat) %in% tss.keep
  return(tss.mat[jkeep, ])
})

print("Dimensions of TSS after keeping one TSS for each gene, defined by highest expression in H3K4me3")
lapply(tss.mats.singlecell, length)




# Get common rows ---------------------------------------------------------

lapply(tss.exprs.lst.unfilt, length)

tss.all <- lapply(tss.exprs.lst, function(exprs.lst){
  names(exprs.lst)
}) %>%
  unlist() %>%
  unique()

tss.common <- lapply(tss.exprs.lst, function(exprs.lst){
  names(exprs.lst) 
}) %>%
  Reduce(f = intersect, .)

# get ensembl names ? 
genes.common <- sapply(tss.common, function(x) strsplit(x, ";")[[1]][[2]])
ens.common <- Gene2Ensembl.ZF(genes.common, return.original = TRUE, species = "drerio")

# create tss, genes, ens dat
genes.annot <- data.frame(bin = tss.common, gene = genes.common, ens = ens.common, stringsAsFactors = FALSE)

# Load UMAP and GLM annots ------------------------------------------------


dat.annot.lst.WKM <- lapply(jmarks, function(jmark){
  print(jmark)
  # filter by previously defined celltypes? 
  inf.annot.louv <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/LDA_downstream/LDA_downstream_ZF.2020-04-23.imputevarfilt.lessstringent/ZF_LDA_output.", jmark, ".keepn_150.final.ClusterTables.txt")
  assertthat::assert_that(file.exists(inf.annot.louv))
  
  
  inf.annot.glmpca <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/celltyping/from_louvain/cell_to_cluster_table.", jmark, ".2020-04-13.txt")
  assertthat::assert_that(file.exists(inf.annot.glmpca))
  
  annot.louv <- fread(inf.annot.louv)
  annot.louv$clusterplate <- paste(annot.louv$cluster, annot.louv$plate, "_")
  
  annot.glmpca <- fread(inf.annot.glmpca)
  annot.glmpca.filt <- subset(annot.glmpca, cell %in% annot.louv$cell) %>%
    rowwise() %>%
    mutate(clusterplate = paste(cluster, plate, sep = "_")) %>%
    mutate(cluster.new = ifelse(cluster %in% c("lymph1", "lymph2"), "lymph", cluster)) %>%   # rename lymph1 and lymph2 into lymph
    ungroup() %>%
    filter(cluster.new != "Unknown") %>%
    left_join(., subset(annot.louv, select = c(cell, var.imputed)))
  # rename clusters
  annot.glmpca.filt$cluster.new <- sapply(annot.glmpca.filt$cluster.new, function(jclst) AssignHash(jclst, wkm.rename, null.fill = jclst))
  print("annot glmpca filt")
  print(annot.glmpca.filt)
  return(annot.glmpca.filt)
})


# filter so only "granu", "HSPCs", "lymph", and "eryth"

dat.annot.lst.WKM <- lapply(dat.annot.lst.WKM, function(jannot){
  jannot$cluster.new <- sapply(jannot$cluster.new, function(x) gsub("monocyte", "granu", x = x))
  jannot$cluster.new <- sapply(jannot$cluster.new, function(x) gsub("HSC", "HSPCs", x = x))
  jannot$cluster.new <- factor(jannot$cluster.new, levels = c("HSPCs", "granu", "lymph", "eryth"))
  jannot <- subset(jannot, !is.na(cluster.new))
  return(jannot)
})


# Set up cells to keep ----------------------------------------------------


cells.keep.lst <- lapply(dat.annot.lst.WKM, function(x){
  return(x$cell)
})



# filter for good cells and good bins
tss.mats.sc.filt.zf <- lapply(jmarks, function(jmark){
  jmat <- tss.mats.singlecell[[jmark]]
  cells.keep <- cells.keep.lst[[jmark]]
  jmat.filt <- jmat[tss.common, cells.keep]
  return(jmat.filt)
})


ncuts.cells <- lapply(tss.mats.singlecell, function(jmat){
  ncuts.dat <- data.frame(cell = colnames(jmat), ncuts.total = colSums(jmat), stringsAsFactors = FALSE)
})

g2e <- hash::hash(genes.common, ens.common)

genesets <- list()

inf.zf.de <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/pseudobulk_analysis.get_DE_genes_from_pbulk_scrnaseq/de_genes_sorted_and_giladi.WithHouseKeepAndNotExpressed.FixExprsOther.exprsmax_5.logfcmin_2.logfcmax_1.RData"
load(inf.zf.de, v=T)

genesets <- de.ens.zf.stringent



jnames <- names(genesets)
names(jnames) <- jnames

genesets.dat <- lapply(jnames, function(jname){
  jgenes <- genesets[[jname]]
  data.frame(geneset = jname, ens = genesets[[jname]], stringsAsFactors = FALSE)
}) %>%
  bind_rows()


print(Sys.time() - jstart)

# save objects
save(tss.mats.sc.filt.zf, ncuts.cells, genesets, genesets.dat, dat.annot.lst.WKM, file = outrdata)


