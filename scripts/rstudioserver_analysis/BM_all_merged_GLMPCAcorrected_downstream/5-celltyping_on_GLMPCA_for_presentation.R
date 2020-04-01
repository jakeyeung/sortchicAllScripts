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
library(forcats)

# jmark <- "H3K4me3"
# jmark <- "H3K4me1"
jmark <- "H3K27me3"
# jmark <- "H3K9me3"

# jcluster <- "HSCs-Hlf_topic26"
# jcluster <- "HSCs-Hlf_topic7"
jcluster <- "HSCs-Tead1-_topic9"
# jcluster <- "HSCs_topic22"
jtopic <- strsplit(jcluster, split = "_")[[1]][[2]]

keeptop <- 20
jprefix <- "/home/jyeung/hub_oudenaarden"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime"
outpdf <- file.path(outdir, paste0("GLMPCA_with_topics.", Sys.Date(), ".", jmark, ".keeptop_", keeptop, ".pdf"))


# Functions ---------------------------------------------------------------

SumAcrossClusters <- function(count.mat, cnames.keep.lst){
  count.mat <- as.matrix(count.mat)
  count.vecs <- lapply(cnames.keep.lst, function(cnames.keep){
    cnames.keep.i <- which(colnames(count.mat) %in% cnames.keep)
    assertthat::assert_that(length(cnames.keep.i) > 0)
    rowSums(count.mat[, cnames.keep.i])
  })
  return(count.vecs)
}



# Load data  --------------------------------------------------------------


print(jmark)

jexperi <- "AllMerged"

mergesize <- "1000"
nbins <- "1000"
jcovar.cname <- "ncuts.var.log2.CenteredAndScaled"
jpenalty <- 1

inf.glm <- file.path(jprefix, paste0("jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2/PZ_", jmark, ".", jexperi, ".KeepBestPlates2.GLMPCA_var_correction.mergebinsize_", mergesize, ".binskeep_", nbins, ".covar_", jcovar.cname, ".penalty_", jpenalty, ".winsorize_TRUE.2020-02-11.RData"))
assertthat::assert_that(file.exists(inf.glm))
load(inf.glm, v=T)


inf.glmpca.annot <- file.path(jprefix, paste0("jyeung/data/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData"))
assertthat::assert_that(file.exists(inf.glmpca.annot))
load(inf.glmpca.annot, v=T)

jinf.tss <- file.path(jprefix, paste0("jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed"))
assertthat::assert_that(file.exists(jinf.tss))

inf.annot <- file.path(jprefix, paste0("jyeung/data/scChiC/public_data/Giladi_et_al_2018/giladi_pseudobulk_datsum_and_ncells.DESeq_and_QuantNorm.RData"))
load(inf.annot, v=T)


inf.lda <- paste0(jprefix, "/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj")
assertthat::assert_that(file.exists(inf.lda))
load(inf.lda, v=T)

tm.result <- posterior(out.lda)
tm.result <- AddTopicToTmResult(tm.result, jsep = "")

topics.mat <- tm.result$topics
# colnames(topics.mat) <- paste("topic", colnames(topics.mat), sep = "")
topics.mat <- data.frame(cell = rownames(topics.mat), topics.mat)


dat.sum.long <- data.table::melt(dat.sum.norm.quantnorm)
colnames(dat.sum.long) <- c("gene", "celltype", "exprs")

dat.sum.long <- dat.sum.long %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE)) %>%
  filter(!is.nan(zscore))

annot.out <- AnnotateCoordsFromList(coords.vec = rownames(glm.out$loadings), inf.tss = jinf.tss)

# Load annotations  -------------------------------------------------------

# findloading for Neutrophils for examle
unique(dat.umap.glm.fillNAs$cluster)



jclust.vec <- as.numeric(dat.umap.glm.fillNAs$cluster == jcluster)
jclust.vec[is.na(jclust.vec)] <- 0

assertthat::assert_that(all(rownames(glm.out$factors) == dat.umap.glm.fillNAs$cell))
glm.out$factors

jcors <- apply(glm.out$factors, MARGIN = 2, FUN = function(jcol){
  cor(jclust.vec, jcol, method = "spearman")
})

dat.umap.glm.merged <- left_join(dat.umap.glm.fillNAs, data.frame(cell = rownames(glm.out$factors), glm.out$factors, stringsAsFactors = FALSE)) %>%
  rowwise() %>%
  mutate(cond = GetCondFromSamp(cell, mark = jmark))

pdf(outpdf, width = 1440/72, height = 815/72, useDingbats = FALSE)


(jdim <- names(jcors)[which.max(abs(jcors))])
print(jcors)
barplot(jcors)

# plot top 
ggplot(dat.umap.glm.merged, aes_string(x = "umap1", y = "umap2", color = jdim)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~cond)

sum(apply(glm.out$factors[, c(1,2)], MARGIN = 1, prod))

sum(apply(glm.out$loadings[, c(1,2)], MARGIN = 1, prod))

# SVD an example
pca.out <- prcomp(glm.out$factors, center = TRUE, scale. = FALSE)
sum(apply(pca.out$rotation[, c(1, 2)], MARGIN = 1, prod))


# Plot LDA loading onto glmpca umap  --------------------------------------



dat.umap.glm.merged.lda <- left_join(dat.umap.glm.fillNAs, topics.mat)
ggplot(dat.umap.glm.merged.lda, aes_string(x = "umap1", y = "umap2", color = jtopic)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c()

# take top 200 bins in LDA and use it to define signature for GLMPCA
dat.impute <- t(as.matrix(glm.out$factors) %*% as.matrix(t(glm.out$loadings)))

# keep topn
bins.keep <- names(sort(tm.result$terms[jtopic, ], decreasing = TRUE)[1:keeptop])

exprs.vec <- colMeans(dat.impute[bins.keep, ])

exprs.dat <- data.frame(cell = colnames(dat.impute), exprs = exprs.vec, stringsAsFactors = FALSE)
dat.umap.glm.impute <- left_join(dat.umap.glm.fillNAs, exprs.dat)
# if we average across many regions maybe we can get something interesting? 

m <- ggplot(dat.umap.glm.impute, aes(x = umap1, y = umap2, color = exprs)) + geom_point() + scale_color_viridis_c(0) + 
  scale_color_viridis_c(direction = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_blank()) + 
  ggtitle(paste(jcluster)) + 
  facet_wrap(~cond)
print(m)

m.merge <- ggplot(dat.umap.glm.impute, aes(x = umap1, y = umap2, color = exprs)) + geom_point() + scale_color_viridis_c(0) + 
  scale_color_viridis_c(direction = 1) + 
  theme_minimal() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text = element_blank()) + 
  ggtitle(paste(jcluster))
print(m.merge)
  
dev.off()

