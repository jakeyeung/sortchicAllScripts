# Jake Yeung
# Date of Creation: 2020-02-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged_GLMPCAcorrected_downstream/1-summarize_GLMPCA_UMAPs.R
# Define celltypes, plot marker genes 


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)
library(scchicFuncs)

library(JFuncs)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(hash)
library(igraph)
library(umap)


# Constants ---------------------------------------------------------------

jmark <- "H3K4me1"
jexperi <- "Unenriched"

nbins <- "500"

jinf.tss <- "/home/jyeung/hpc/databases/gene_tss/gene_tss_winsize.50000.bed"
assertthat::assert_that(file.exists(jinf.tss))

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# Load LDA output ---------------------------------------------------------


inf.lda <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-11.", jexperi, ".K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-11.", jexperi, ".K-30.Robj")
assertthat::assert_that(file.exists(inf.lda))
load(inf.lda, v=T)

tm.result <- posterior(out.lda)

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

annots.out <- AnnotateBins2(tm.result$terms, top.thres = 0.995, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db")
annots.out$terms.annot$gene <- sapply(annots.out$terms.annot$termgene, function(x) ifelse(is.na(x), NA, strsplit(x, ";")[[1]][[2]]))

termgene.dat <- annots.out$terms.annot %>%
  group_by(term) %>%
  filter(rnk == min(rnk))

dat.umap.lda <- DoUmapAndLouvain(topics.mat = posterior(out.lda)$topics, jsettings = jsettings)

# Load GLMPCA correctoin --------------------------------------------------


jcovar.cname <- "ncuts.var"; inf.glm <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2/PZ_", jmark, ".", jexperi, ".KeepBestPlates2.GLMPCA_var_correction.mergebinsize_1000.binskeep_", nbins, ".devmode.2020-02-11.RData")


jcovar.cname <- "cell.var.within.sum.norm.CenteredAndScaled"
jcovar.cname <- "ncuts.var.CenteredAndScaled"
# inf.glm <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.KeepBestPlates2/PZ_", jmark, ".", jexperi, ".KeepBestPlates2.GLMPCA_var_correction.mergebinsize_1000.binskeep_250.covar_", jcovar.cname, ".devmode.2020-02-11.RData")
assertthat::assert_that(file.exists(inf.glm))
load(inf.glm, v=T)

tm.result.glm <- list(topics = as.matrix(glm.out$factors), terms = t(as.matrix(glm.out$loadings)))
dat.impute.glm <- t(tm.result.glm$topics %*% tm.result.glm$terms)

# Pload umap before and after  ---------------------------------------------


cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

dat.umap.glm <- DoUmapAndLouvain(topics.mat = glm.out$factors, jsettings = jsettings)

m.lda <- ggplot(dat.umap.lda, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = cbPalette) + ggtitle("From LDA")
m.glm <- ggplot(dat.umap.glm, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = cbPalette) + ggtitle("GLMPCA")

multiplot(m.lda, m.glm, cols = 2)

# check batch correction loadings?
batch.cors <- data.frame(term = rownames(glm.out$coefX), 
                         batch.cor = glm.out$coefX[[jcovar.cname]], stringsAsFactors = FALSE) %>%
  left_join(., termgene.dat)

plot(density(batch.cors$batch.cor))

# show a term with large loading
jterm <- subset(batch.cors, batch.cor == min(batch.cor))$term[[1]]
jsub <- subset(batch.cors, gene == "S100a8")
jterm <- jsub$term[[1]]

# plot this gene in original LDA
corrected.vec <- data.frame(cell = colnames(dat.impute.glm), exprs.corrected = dat.impute.glm[jterm, ], stringsAsFactors = FALSE)
impute.vec <- data.frame(cell = colnames(dat.impute.log), exprs = dat.impute.log[jterm, ], stringsAsFactors = FALSE)

jdat <- left_join(dat.umap.lda, impute.vec) %>%
  left_join(., subset(dat.merge2, select = c(cell, cell.var.within.sum.norm, ncuts.var, ncuts))) %>%
  left_join(., corrected.vec)

jdat$correction <- jdat$ncuts.var * jsub$batch.cor

m.check <- ggplot(jdat, aes(x = umap1, y = umap2, color = exprs)) + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_viridis_c() + 
  ggtitle("From LDA")
print(m.check)

m.correction <- ggplot(jdat, aes(x = umap1, y = umap2, color = correction)) + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + 
  ggtitle("Correction vec")
print(m.correction)

m.corrected <- ggplot(jdat, aes(x = umap1, y = umap2, color = exprs.corrected)) + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.position = "bottom") + 
  scale_color_viridis_c() + 
  ggtitle("correccted expression from glmpca")
print(m.corrected)

multiplot(m.check, m.corrected, cols = 2)

m.check.lda <- ggplot(jdat, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = cbPalette) + 
  ggtitle("From LDA")
print(m.check.lda)

m.before <- ggplot(jdat, aes(x = exprs, y = cell.var.within.sum.norm, color = louvain)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)
m.after <- ggplot(jdat, aes(x = exprs.corrected, y = cell.var.within.sum.norm, color = louvain)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_manual(values = cbPalette)
multiplot(m.before, m.after)

ggplot(jdat, aes(x = exprs, y = exprs.corrected, color = louvain)) + geom_point() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# ggplot(jdat, aes(x = exprs, y = ncuts.var)) + geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# plot corrected
 
# Label celltypes using topics --------------------------------------------






