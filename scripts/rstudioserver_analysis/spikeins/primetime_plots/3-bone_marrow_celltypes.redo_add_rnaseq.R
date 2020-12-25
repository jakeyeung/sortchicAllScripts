# Jake Yeung
# Date of Creation: 2020-12-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/3-bone_marrow_celltypes.redo_add_rnaseq.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(topicmodels)

library(hash)
library(igraph)
library(umap)


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
jmarks.act <- c("H3K4me1", "H3K4me3"); names(jmarks.act) <- jmarks.act
jmarks.rep <- c("H3K27me3" = "H3K27me3")

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")


# Load glmpcas  -----------------------------------------------------------


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Load annots  ------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"


inf.meta.k4me1 <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.H3K4me1.txt"))
inf.meta.k4me3 <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-04.more_celltypes.metadta.H3K4me3.txt"))
inf.meta.k27me3 <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.H3K27me3_techrepmerged/BM_rep2_rep3reseq_H3K27me3.2020-12-10.txt"))

dat.k4me1 <- fread(inf.meta.k4me1)
dat.k4me3 <- fread(inf.meta.k4me3)
dat.k27me3 <- fread(inf.meta.k27me3) %>%
  dplyr::rename(cluster.more = cluster,
                cluster = cluster.fewer) 

ggplot(dat.k4me1, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.k4me3, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.k27me3, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())




# Load gene sets ----------------------------------------------------------

inf.gsets <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/pdfs_all/genesets_on_umap/geneset_on_umap.binskeep_0.niter_500.2020-12-09.from_LDA_topics.condensed.heatmap.famousgenes.keepn_400.refmark_H3K4me3.2020-12-09.txt")

dat.genes <- fread(inf.gsets)

# get coordinate which will be used to rename K27me3
dat.genes$coord <- sapply(dat.genes$gene, function(x) paste("chr", strsplit(x, ";")[[1]][[1]], sep = ""))

dat.genes <- as.data.frame(dat.genes)

dat.genes$symbol <- sapply(dat.genes$gene, function(g) strsplit(g, split = "\\.")[[1]][[4]])

dat.genes <- dat.genes %>%
  group_by(jset) %>%
  mutate(rnk = seq(length(gene))) 

dat.genes.filt <- subset(dat.genes, rnk <= 150)


# Add colors to annotation ------------------------------------------------

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
jsub.annot.new <- dat.k4me1 %>% arrange(cluster, jrep)
clstrs.uniq.new <- sort(unique(dat.genes$jset))
batch.uniq.new <- unique(jsub.annot.new$jrep)
colors.uniq.new <- cbPalette[1:length(clstrs.uniq.new)]
colors.uniq.batch.new <- cbPalette[10:(9 + length(unique(jsub.annot.new$jrep)))]
cols.hash.new <- hash::hash(clstrs.uniq.new, colors.uniq.new)
cols.hash.batch.new <- hash::hash(batch.uniq.new, colors.uniq.batch.new)
colsvec.new <- sapply(jsub.annot.new$cluster2, function(x) AssignHash(x, cols.hash.new, null.fill = NA))
colsvec.batch.new <- sapply(jsub.annot.new$jrep, function(x) AssignHash(x, cols.hash.batch.new, null.fill = NA))
colsvec.row.new <- sapply(dat.genes$jset, function(x) AssignHash(x, cols.hash.new, null.fill = NA))

dat.meta.new.lst <- list("H3K4me1" = dat.k4me1, "H3K4me3" = dat.k4me3, "H3K27me3" = dat.k27me3)
dat.meta.new.lst <- lapply(dat.meta.new.lst, function(jdat){
  jdat$colorcode <- sapply(jdat$cluster, function(x) AssignHash(x, cols.hash.new, null.fill = NULL))
  jdat <- jdat %>%
    arrange(cluster, jrep)
  return(jdat)
})

# color genes
dat.genes$colorcode <- sapply(dat.genes$jset, function(x) AssignHash(x, cols.hash.new, null.fill = NULL))


genes2col <- hash::hash(dat.genes$gene, dat.genes$colorcode)

jtest <- "H3K27me3"
ggplot(dat.meta.new.lst[[jtest]], aes(x = umap1, y = umap2, color = colorcode)) + 
  geom_point() + 
  theme_bw() + 
  scale_color_identity() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load raw data  ----------------------------------------------------------

infs.raw.act.lst <- lapply(jmarks.act, function(jmark){
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.from_TSS.dist_10000.RemoveDupRows/lda_outputs.count_mat_from_TSS.", jmark, ".dist_10000.K-30.binarize.FALSE/ldaOut.count_mat_from_TSS.", jmark, ".dist_10000.K-30.Robj"))
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

infs.raw.rep.lst <- lapply(jmarks.rep, function(jmark){
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_rep2rep3reseq_TSS/lda_outputs.PZ-BM-rep3-", jmark, "-rep2rep3reseq.TSS.varfilt.K-30.binarize.FALSE/ldaOut.PZ-BM-rep3-", jmark, "-rep2rep3reseq.TSS.varfilt.K-30.Robj"))
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

infs.raw.lst <- c(infs.raw.act.lst, infs.raw.rep.lst)

out.lda.lst <- lapply(infs.raw.lst, function(x){
  load(x, v=T)
  return(list(count.mat = count.mat, out.lda = out.lda))
})

tm.result.lst <- lapply(out.lda.lst, function(x){
  tm.result <- posterior(x$out.lda)
  tm.result <- AddTopicToTmResult(tm.result, jsep = "")
  return(tm.result)
})


# Load raw K27me3  --------------------------------------------------------




# Rename rownames for K27me3  ---------------------------------------------

rnames.act <- lapply(jmarks.act, function(jmark){
  rnames <- rownames(out.lda.lst[[jmark]]$count.mat)
}) %>%
  unlist() %>%
  unique()  %>%
  gtools::mixedsort()

coords.act <- sapply(rnames.act, function(x) paste("chr", strsplit(x, split = ";")[[1]][[1]], sep = ""))

coord2rname <- hash::hash(coords.act, rnames.act)


# Make heatmaps with raw  ----------------------------------------------------------

count.mat.norm.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  count.mat.norm <- sweep(out.lda.lst[[jmark]]$count.mat, MARGIN = 2, STATS = colSums(out.lda.lst[[jmark]]$count.mat), FUN = "/")
  # if h3k27me3, rename rownames
  if (jmark == "H3K27me3"){
    print("Renaming rownames")
    rnames.before <- rownames(count.mat.norm) 
    rnames.after <- sapply(rnames.before, function(coord) AssignHash(x = coord, jhash = coord2rname, null.fill = NA))
    rnames.filt <- rnames.after[!is.na(rnames.after)]
    count.mat.norm <- count.mat.norm[names(rnames.filt), ]
    # rename
    rownames(count.mat.norm) <- rnames.filt
  }
  return(count.mat.norm)
})

count.mat.raw.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  count.mat <- out.lda.lst[[jmark]]$count.mat
  # if h3k27me3, rename rownames
  if (jmark == "H3K27me3"){
    print("Renaming rownames")
    rnames.before <- rownames(count.mat) 
    rnames.after <- sapply(rnames.before, function(coord) AssignHash(x = coord, jhash = coord2rname, null.fill = NA))
    rnames.filt <- rnames.after[!is.na(rnames.after)]
    count.mat <- count.mat[names(rnames.filt), ]
    rownames(count.mat) <- rnames.filt
  }
  return(count.mat)
})

# rename for tm.result k27me3

colnames(tm.result.lst$H3K27me3$terms) <- sapply(colnames(tm.result.lst$H3K27me3$terms), function(coord) AssignHash(x = coord, jhash = coord2rname, null.fill = coord))


# Make heatmap ------------------------------------------------------------



jtest <- "H3K4me3"

jtest <- "H3K4me1"

count.mat.tmp <- count.mat.norm.lst[[jtest]]

gkeep <- dat.genes$gene
rkeep <- sapply(dat.genes$gene, function(g) paste("chr", strsplit(g, ";")[[1]][[1]], sep = ""))

# not.in.rows <- rkeep[which(!rkeep %in% rownames(count.mat.tmp))]
# in.rows <- rkeep[which(rkeep %in% rownames(count.mat.tmp))]

ckeep <- dat.meta.new.lst[[jtest]]$cell
rowcols <- dat.genes$colorcode
colcols <- dat.meta.new.lst[[jtest]]$colorcode

# gkeep.i <- rownames(count.mat.tmp) %in% gkeep
# ckeep.i <- colnames(count.mat.tmp) %in% ckeep
# # count.mat.tmp <- count.mat.norm.lst[[jtest]][rkeep, ckeep]
# count.mat.tmp <- count.mat.norm.lst[[jtest]][gkeep.i, ckeep.i]

# count.mat.tmp.filt <- count.mat.tmp[rkeep, ckeep]
count.mat.tmp.filt <- count.mat.tmp[gkeep, ckeep]

print(dim(count.mat.tmp.filt))
# assertthat::assert_that(length(gkeep.i) > 0)
# assertthat::assert_that(length(ckeep.i) > 0)

heatmap3::heatmap3(as.matrix(BinarizeMatrix(count.mat.tmp.filt)), Rowv = NA, Colv = NA, ColSideColors = colcols, scale = "row", RowSideColors = rowcols, revC = TRUE, main = paste(jtest))



# Compare gene sets with Giladi  ------------------------------------------

m2c <- MarkerToCelltype()
inf.giladi <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Giladi_et_al_2018/giladi_pseudobulk_datsum_and_ncells.DESeq_and_QuantNorm.RData")
load(inf.giladi, v=T)

dat.sum.long <- data.table::melt(dat.sum.norm.quantnorm)
colnames(dat.sum.long) <- c("gene", "celltype", "exprs")

dat.sum.long <- dat.sum.long %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE)) %>%
  filter(!is.nan(zscore)) 


dat.sum.mat <- reshape2::dcast(dat.sum.long, formula = "gene ~ celltype", value.var = "exprs") %>%
# dat.sum.mat <- reshape2::dcast(dat.sum.long, formula = "gene ~ celltype", value.var = "zscore") %>%
  as.data.frame()

rownames(dat.sum.mat) <- dat.sum.mat$gene
dat.sum.mat$gene <- NULL

colnames(dat.sum.mat) <- sapply(colnames(dat.sum.mat), function(x) m2c[[x]])
dat.sum.mat <- dat.sum.mat[, order(colnames(dat.sum.mat))]

# genes.grep <- paste(unique(dat.genes.filt$symbol), collapse = "|")

gvec <- dat.genes.filt$symbol
names(gvec) <- dat.genes.filt$gene





gvec <- subset(dat.genes, rnk <= 300)$symbol
names(gvec) <- subset(dat.genes, rnk <= 300)$gene

gvec.uniq <- gvec[!duplicated(gvec)]

genes.grep.i.lst <- lapply(gvec.uniq, function(x) grep(x, rownames(dat.sum.mat)))

# remove no matches
genes.grep.i <- genes.grep.i.lst[unlist(lapply(genes.grep.i.lst, function(x) length(x) > 0))]

# remove dupes
genes.grep.i.dedup <- sapply(genes.grep.i, function(x) x[[1]])

dat.sum.mat.filt <- dat.sum.mat[unlist(genes.grep.i.dedup), ]

colors.genes <- sapply(names(genes.grep.i.dedup), function(x) AssignHash(x, genes2col, null.fill = x))

heatmap3::heatmap3(dat.sum.mat.filt, Rowv = NA, Colv = NA, scale = "row", RowSideColors = colors.genes, revC = TRUE, main = paste("Giladi"))

# heatmap3::heatmap3(dat.sum.mat.filt, Rowv = NA, Colv = NA, ColSideColors = NULL, scale = "row", RowSideColors = colors.genes, revC = TRUE, main = paste("Giladi"))



# Write outputs share with others -----------------------------------------

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/with_imputed"
dir.create(outdir)

# write metadat
for (jmark in jmarks){
  fnametmp <- paste0("metadata_umap_celltype_cuts.", jmark, ".txt")
  outtmp <- file.path(outdir, fnametmp)
  fwrite(dat.meta.new.lst[[jmark]], file = outtmp)
}

# write imputed
outtmp.imputed <- file.path(outdir, "imputed_TSS_three_marks.rds")
saveRDS(tm.result.lst, outtmp.imputed)

# write Giladi object
outtmp.giladi <- file.path(outdir, "giladi_pseudobulk_exprs_data.rds")
saveRDS(dat.sum.long, outtmp.giladi)

# write raw counts at TSS
outtmp.raw <- file.path(outdir, "raw_cuts_TSS_three_marks.rds")
saveRDS(count.mat.raw.lst, outtmp.raw)

# write gene sets
fwrite(dat.genes, file = file.path(outdir, "celltype_specific_genes_defined_by_K4me3_TSS.txt"))


# Make heatmaps for Alexander ---------------------------------------------

jtest <- "H3K4me3"
jtest <- "H3K4me1"

outpdf <- file.path(outdir, "heatmap_three_marks_gene_set.pdf")
pdf(outpdf, useDingbats = FALSE)
for (jtest in jmarks){
  print(jtest)
  count.mat.tmp <- count.mat.norm.lst[[jtest]]
  
  gkeep <- dat.genes$gene
  rkeep <- sapply(dat.genes$gene, function(g) paste("chr", strsplit(g, ";")[[1]][[1]], sep = ""))
  
  # not.in.rows <- rkeep[which(!rkeep %in% rownames(count.mat.tmp))]
  # in.rows <- rkeep[which(rkeep %in% rownames(count.mat.tmp))]
  
  ckeep <- dat.meta.new.lst[[jtest]]$cell
  rowcols <- dat.genes$colorcode
  colcols <- dat.meta.new.lst[[jtest]]$colorcode
  
  count.mat.tmp.filt <- count.mat.tmp[gkeep, ckeep]
  
  print(dim(count.mat.tmp.filt))
  print(count.mat.tmp.filt[1:5, 1:5])
  
  heatmap3::heatmap3(as.matrix(BinarizeMatrix(count.mat.tmp.filt)), Rowv = NA, Colv = NA, ColSideColors = colcols, scale = "row", RowSideColors = rowcols, revC = TRUE, main = paste(jtest, "cuts"))
  
  # sum across cells
  genes.lst <- split(x = dat.genes$gene, f = dat.genes$jset)
  clstrs.lst <- split(x = dat.meta.new.lst[[jtest]]$cell, f = dat.meta.new.lst[[jtest]]$cluster)
  
  pseudogene.mat <- SumAcrossClusters(t(count.mat.tmp.filt), cnames.keep.lst = genes.lst)
  pseudogene.mat <- do.call(rbind, pseudogene.mat)
  
  pseudocell.mat <- SumAcrossClusters(count.mat.tmp.filt, cnames.keep.lst = clstrs.lst)
  pseudocell.mat <- do.call(cbind, pseudocell.mat)
  
  
  # sum across genes
  heatmap3::heatmap3(log2(pseudogene.mat + 0.01), Rowv = NA, Colv = NA, ColSideColors = colcols, scale = "row", RowSideColors = sapply(rownames(pseudogene.mat), AssignHash, cols.hash.new), revC = TRUE, main = paste(jtest, "pseudogene"))
  heatmap3::heatmap3(log2(pseudocell.mat + 0.01), Rowv = NA, Colv = NA, ColSideColors = sapply(colnames(pseudocell.mat), AssignHash, cols.hash.new), scale = "row", RowSideColors = rowcols, revC = TRUE, main = paste(jtest, "pseudobulk"))
}
dev.off()



