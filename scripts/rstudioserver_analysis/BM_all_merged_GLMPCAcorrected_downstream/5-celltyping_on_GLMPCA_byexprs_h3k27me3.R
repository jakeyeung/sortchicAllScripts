# Jake Yeung
# Date of Creation: 2020-02-19
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged_GLMPCAcorrected_downstream/5-celltyping_on_GLMPCA.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(topicmodels)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(scchicFuncs)


# Functions ---------------------------------------------------------------


# Load data  --------------------------------------------------------------

jmark <- "H3K27me3"
print(jmark)

jexperi <- "AllMerged"

mergesize <- "1000"
nbins <- "1000"
jcovar.cname <- "ncuts.var.log2.CenteredAndScaled"
jpenalty <- 1

inf.glm <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.KeepBestPlates2/PZ_", jmark, ".", jexperi, ".KeepBestPlates2.GLMPCA_var_correction.mergebinsize_", mergesize, ".binskeep_", nbins, ".covar_", jcovar.cname, ".penalty_", jpenalty, ".winsorize_TRUE.2020-02-11.RData")
assertthat::assert_that(file.exists(inf.glm))
load(inf.glm, v=T)


inf.glmpca.annot <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
assertthat::assert_that(file.exists(inf.glmpca.annot))
load(inf.glmpca.annot, v=T)

jinf.tss <- "/home/jyeung/hpc/databases/gene_tss/gene_tss_winsize.50000.bed"
assertthat::assert_that(file.exists(jinf.tss))

inf.annot <- "/home/jyeung/hpc/scChiC/public_data/Giladi_et_al_2018/giladi_pseudobulk_datsum_and_ncells.DESeq_and_QuantNorm.RData"
load(inf.annot, v=T)

dat.sum.long <- data.table::melt(dat.sum.norm.quantnorm)
colnames(dat.sum.long) <- c("gene", "celltype", "exprs")

dat.sum.long <- dat.sum.long %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE)) %>%
  filter(!is.nan(zscore))


# Compare loadings with public data ---------------------------------------

# get common rows then plot correlations?

inf.tss <- "/home/jyeung/hpc/databases/gene_tss/gene_tss_winsize.50000.bed"
annot.out <- AnnotateCoordsFromList(coords.vec = rownames(glm.out$loadings), inf.tss = inf.tss)

# x <- glm.out$loadings[, jdim]
# rename region_coord to symbol?

rnames.old <- rownames(glm.out$loadings)
coord2gene <- hash::hash(annot.out$regions.annotated$region_coord, annot.out$regions.annotated$SYMBOL)

glm.out$loadings[, jdim]

loadings.renamed <- data.frame(coord = rownames(glm.out$loadings), 
                               gene = sapply(rownames(glm.out$loadings), AssignHash, jhash = coord2gene),
                               glm.out$loadings, stringsAsFactors = FALSE)

# add gene expression (zscore) to the loadings?
loadings.sub <- loadings.renamed[, c("coord", "gene", jdim)]
# add zscore across cells
dat.sum.zscore <- reshape2::dcast(data = dat.sum.long, formula = gene ~ celltype, value.var = "zscore")
# loadings.sub <- left_join(loadings.sub, dat.sum.zscore)

#  get common genesD
genes.common <- intersect(dat.sum.zscore$gene, loadings.sub$gene)

dat.sum.zscore.merge <- inner_join(dat.sum.zscore, loadings.sub) %>%
  filter(!duplicated(gene))

ctypes <- colnames(dat.sum.norm)


# Plot imputed expression of marker genes? --------------------------------


dat.impute <- t(as.matrix(glm.out$factors) %*% as.matrix(t(glm.out$loadings)))

# jgene <- "S100a7"
jgene <- "Bach2"
jgene <- "Siglech"
jgene <- "Irf4"

jgene <- "Hlf"

jgene <- "Ltf"
(jsub <- subset(annot.out$regions.annotated, grepl(jgene, SYMBOL)))
jterm <- jsub$region_coord[[1]]

exprs.dat <- data.frame(cell = colnames(dat.impute), exprs = dat.impute[jterm, ], stringsAsFactors = FALSE)

dat.umap.glm.impute <- left_join(dat.umap.glm.fillNAs, exprs.dat)

# if we average across many regions maybe we can get something interesting? 
ggplot(dat.umap.glm.impute, aes(x = umap1, y = umap2, color = exprs)) + geom_point() + scale_color_viridis_c(0) + 
  scale_color_viridis_c(direction = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste(jgene, jterm))



# Use topic loadings from K27me3 to plot GLMPCA outputs --------------------

# load active mark

jmark.act <- "H3K4me1"
jmark.repress <- "H3K27me3"
inf.act <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jmark.act, "_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jmark.act, "_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj")
assertthat::assert_that(file.exists(inf.act))

inf.repress <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jmark.repress, "_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jmark.repress, "_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj")
assertthat::assert_that(file.exists(inf.repress))

load(inf.act, v=T)
count.mat.act <- count.mat
out.lda.act <- out.lda

tm.result.act <- posterior(out.lda)
rownames(tm.result.act$terms) <- paste("topic", rownames(tm.result.act$terms), sep = "")
colnames(tm.result.act$topics) <- paste("topic", colnames(tm.result.act$topics), sep = "")


load(inf.repress, v=T)
count.mat.repress <- count.mat
out.lda.repress <- out.lda

tm.result.repress <- posterior(out.lda)
rownames(tm.result.repress$terms) <- paste("topic", rownames(tm.result.repress$terms), sep = "")
colnames(tm.result.repress$topics) <- paste("topic", colnames(tm.result.repress$topics), sep = "")


# Take top genes by zscore and plot in UMAP  ------------------------------



print(unique(as.character(dat.sum.long$celltype)))

jctype <- "core"
jctype <- "Ccl5"
jctype <- "Fcrla"
jctype <- "Prss34"
jctype <- "Ltf"

pdf("/home/jyeung/hpc/scChiC/from_rstudioserver/pdfs_all/stemcell_analysis/celltype_specific_genes_on_K27me3_signal.pdf", useDingbats = FALSE)
for (jctype in unique(as.character(dat.sum.long$celltype))){
  print(jctype)
  
  jsub.ref <- subset(dat.sum.long, celltype == jctype) %>%
    arrange(desc(zscore)) %>%
    # arrange(zscore) %>%
    ungroup() %>%
    mutate(gene = as.character(gene),
           rnk = rank(-zscore, ties.method = "min"))
  
  jgenes <- jsub.ref$gene[1:500]
  # jgenes <- sample(subset(dat.sum.long, celltype == jctype)$gene, size = 100, replace = FALSE)
 
  # jgenes.all <- annot.out$regions.annotated$SYMBOL
  # (jgenes <- sample(jgenes.all, size = 100, replace = FALSE))
  
  # jgenes <- c("Sez6", "Tead1")
  # coords.keep.hsc <- names(sort(tm.result.repress$terms["topic9", ], decreasing = TRUE)[1:1000])
  # jgenes <- subset(annot.out$regions.annotated, region_coord %in% coords.keep.hsc)$SYMBOL
  
  # plot signal of K4me1 onto umap
  dat.coords.keep.ref <- subset(annot.out$regions.annotated, SYMBOL %in% jgenes) %>%
  # dat.coords.keep.ref <- subset(annot.out$regions.annotated, region_coord %in% coords.keep) %>%
    group_by(SYMBOL) %>%
    filter(abs(distanceToTSS) == min(distanceToTSS)) %>%
    arrange(SYMBOL)
  
  coords.keep <- dat.coords.keep.ref$region_coord
  
  print(coords.keep)
  
  # coords.keep <- sample(colnames(tm.result.repress$term), size = 100, replace = FALSE)
  
  rows.i <- which(rownames(dat.impute) %in% coords.keep)
  print(rows.i)
  # plot in GLMPCA, average across egions 
  jsub <- dat.impute[rows.i, ]
  dat.sub <- data.frame(cell = colnames(dat.impute), exprs = colMeans(jsub), stringsAsFactors = FALSE) %>%
    left_join(., subset(dat.umap.glm.fillNAs, select = c(umap1, umap2, cell, cluster))) %>%
    rowwise() %>%
    mutate(cond = GetCondFromSamp(cell, mark = "H3K27me3"))
  
  # clusters.remove <- c("topic9", "topic26", "topic17", "topic27", "topic30")
  # clusters.remove <- c()
  # clusters.keep <- unique(subset(dat.sub, !grepl(paste(clusters.remove, collapse = "|"), cluster))$cluster)
  
  m <- ggplot(dat.sub, aes(x = umap1, y = umap2, color = exprs)) + geom_point() + scale_color_viridis_c(0) + 
    scale_color_viridis_c(direction = -1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle("Randonly sample 100 genes from genome") + 
    ggtitle(paste(jctype, "avg across", length(coords.keep), "regions"))
  m.split <- ggplot(dat.sub, aes(x = umap1, y = umap2, color = exprs)) + geom_point() + scale_color_viridis_c(0) + 
    scale_color_viridis_c(direction = -1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle("Randonly sample 100 genes from genome") + 
    ggtitle(paste(jctype, "avg across", length(coords.keep), "regions")) + 
    facet_wrap(~cond)
  print(m)
  print(m.split)
}
dev.off()



