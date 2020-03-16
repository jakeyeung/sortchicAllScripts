# Jake Yeung
# Date of Creation: 2020-03-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/stemcell_analysis/check_K27me3_gene_distribution.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)

library(scchicFuncs)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)


# load GLMPCA from bins 
jmark <- "H3K27me3"
jexperi <- "AllMerged"
mergesize <- "1000"
nbins <- "1000"
jcovar.cname <- "ncuts.var.log2.CenteredAndScaled"
jpenalty <- 1
ntopics <- 30

inf.glmpca <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.good_runs/PZ_", jmark, ".AllMerged.KeepBestPlates2.GLMPCA_var_correction.mergebinsize_1000.binskeep_1000.covar_ncuts.var.log2.CenteredAndScaled.penalty_1.winsorize_TRUE.2020-02-11.RData")
inf.lda <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping/GLMPCA_celltyping.", jmark, ".AllMerged.mergesize_1000.nbins_1000.penalty_1.covar_ncuts.var.log2.CenteredAndScaled.RData")
inf.lda.bins <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-11.AllMerged.K-30.Robj")

load(inf.glmpca, v=T)
load(inf.lda, v=T)
load(inf.lda.bins, v=T)

# Annotate bins -----------------------------------------------------------

jinf.tss <- "/home/jyeung/hpc/databases/gene_tss/gene_tss_winsize.50000.bed"
assertthat::assert_that(file.exists(jinf.tss))

tm.result <- posterior(out.lda)

terms.mat <- tm.result$terms
topics.mat <- tm.result$topics

annots.out <- AnnotateBins2(terms.mat = terms.mat, top.thres = 0.995, 
                            inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                            annodb = "org.Mm.eg.db")


# Plot LDA ----------------------------------------------------------------

# get fraction of bins within 50kb of gene vs not
jdist <- 0
prox.bins <- subset(annots.out$regions.annotated, abs(distanceToTSS) <= jdist)$region_coord
dist.bins <- subset(annots.out$regions.annotated, abs(distanceToTSS) > jdist)$region_coord

sample(x = prox.bins, size = 1)
sample(x = dist.bins, size = 1)

prox.fracs <- colSums(count.mat[prox.bins, ]) / colSums(count.mat)
# dist.fracs <- colSums(count.mat[dist.bins, ]) / colSums(count.mat[c(prox.bins, dist.bins), ])

prox.fracs.dat <- data.frame(cell = names(prox.fracs), prox.frac = prox.fracs, stringsAsFactors = FALSE)

dat.umap.annot <- left_join(dat.umap.glm.fillNAs, prox.fracs.dat)

# Show variance corrected fracs -------------------------------------------

dat.impute <- t(as.matrix(glm.out$factors) %*% as.matrix(t(glm.out$loadings)))

rowsd <- apply(dat.impute, 1, sd)
dat.impute.zscore <- sweep(dat.impute, MARGIN = 1, STATS = rowsd, FUN = "/")

jbins <- intersect(rownames(dat.impute), prox.bins)

# prox.fracs.impute <- colSums(dat.impute[jbins, ]) / colSums(dat.impute)
# prox.fracs.impute.dat <- data.frame(cell = names(prox.fracs.impute), prox.frac.impute = prox.fracs.impute, stringsAsFactors = FALSE)

prox.fracs.impute <- colMeans(dat.impute.zscore[jbins, ])
prox.fracs.impute.dat <- data.frame(cell = names(prox.fracs.impute), prox.frac.impute = prox.fracs.impute, stringsAsFactors = FALSE)

dat.umap.annot.impute <- left_join(dat.umap.glm.fillNAs, prox.fracs.impute.dat)





# Plot outputs ------------------------------------------------------------

outpdf <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/pdfs_all/stemcell_analysis/proxfrac_across_cells.", jmark, ".pdf")
pdf(outpdf, useDingbats = FALSE)

ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = prox.frac)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle(jmark, "Fraction of reads in bins that overlap TSS")

ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = prox.frac)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle(jmark, "Fraction of reads in bins that overlap TSS") +
  facet_wrap(~cond)

ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = prox.frac)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + ggtitle(jmark, "Fraction of reads in bins that overlap TSS")

ggplot(dat.umap.annot, aes(x = umap1, y = umap2, color = prox.frac)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + ggtitle(jmark, "Fraction of reads in bins that overlap TSS") + 
  facet_wrap(~cond)

ggplot(dat.umap.annot.impute, aes(x = umap1, y = umap2, color = prox.frac.impute)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle(jmark, "Imputed zscore in bins that overlap TSS")

ggplot(dat.umap.annot.impute, aes(x = umap1, y = umap2, color = prox.frac.impute)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + ggtitle(jmark, "Imputed zscore in bins that overlap TSS") + 
  facet_wrap(~cond)

ggplot(dat.umap.annot.impute, aes(x = umap1, y = umap2, color = prox.frac.impute)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + ggtitle(jmark, "Imputed zscore in bins that overlap TSS")

ggplot(dat.umap.annot.impute, aes(x = umap1, y = umap2, color = prox.frac.impute)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + ggtitle(jmark, "Imputed zscore in bins that overlap TSS") + 
  facet_wrap(~cond)

ggplot(dat.umap.annot.impute, aes(x = prox.frac.impute, fill = cond)) + geom_density(alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste("Imputed zscore", jmark))

# Plot density  -----------------------------------------------------------

ggplot(dat.umap.annot, aes(x = prox.frac, fill = cond)) + geom_density(alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste("Fraction of reads in bins that overlap TSS", jmark))

ggplot(dat.umap.annot.impute, aes(x = prox.frac.impute, fill = cond)) + geom_density(alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste("Imputed zscore", jmark))

dev.off()

