# Jake Yeung
# Date of Creation: 2020-02-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged_GLMPCAcorrected_downstream/2-summarize_GLMPCA_UMAPs.R
# 

rm(list=ls())

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

library(mixtools)


# Constants ---------------------------------------------------------------

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

jthres.frac.cells <- 0.9  # threshold for imputing NAs to celltype 
# (e.g. NA associarted to louvain cluster that are 90% of one celltype)

jmark <- "H3K4me1"
# jmark <- "H3K4me1"
# jmark <- "H3K4me3"
# jmark <- "H3K27me3"
# jmark <- "H3K9me3"
jexperi <- "Unenriched"
# jexperi <- "AllMerged"

mergesize <- "1000"
nbins <- "1000"
# jcovar.cname <- "ncuts.var"; inf.glm <- paste0("/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2/PZ_", jmark, ".", jexperi, ".KeepBestPlates2.GLMPCA_var_correction.mergebinsize_1000.binskeep_", nbins, ".devmode.2020-02-11.RData")
# jcovar.cname <- "cell.var.within.sum.norm.CenteredAndScaled"
jcovar.cname <- "ncuts.var.log2.CenteredAndScaled"
# jcovar.cname <- "ncuts.var.CenteredAndScaled"
jpenalty <- 1
# inf.glm <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.KeepBestPlates2/PZ_", jmark, ".", jexperi, ".KeepBestPlates2.GLMPCA_var_correction.mergebinsize_1000.binskeep_", nbins, ".covar_", jcovar.cname, ".devmode.2020-02-11.RData")
inf.glm <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.KeepBestPlates2/PZ_", jmark, ".", jexperi, ".KeepBestPlates2.GLMPCA_var_correction.mergebinsize_", mergesize, ".binskeep_", nbins, ".covar_", jcovar.cname, ".penalty_", jpenalty, ".winsorize_TRUE.2020-02-11.RData")
assertthat::assert_that(file.exists(inf.glm))

jinf.tss <- "/home/jyeung/hpc/databases/gene_tss/gene_tss_winsize.50000.bed"
assertthat::assert_that(file.exists(jinf.tss))

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

pdfdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping"
pdfname <- paste0("GLMPCA_celltyping.", jmark, ".", jexperi, ".mergesize_", mergesize, ".nbins_", nbins, ".penalty_", jpenalty, ".covar_", jcovar.cname, ".pdf")
pdfout <- file.path(pdfdir, pdfname)

pdf(pdfout, useDingbats = FALSE)
# Load LDA output ---------------------------------------------------------


inf.lda <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-02-11.var_filt.UnenrichedAndAllMerged.KeepBestPlates2/lda_outputs.BM_", jmark, "_varfilt_countmat.2020-02-11.", jexperi, ".K-30.binarize.FALSE/ldaOut.BM_", jmark, "_varfilt_countmat.2020-02-11.", jexperi, ".K-30.Robj")
assertthat::assert_that(file.exists(inf.lda))
load(inf.lda, v=T)

tm.result <- posterior(out.lda)
# label  tpic names
colnames(tm.result$topics) <- paste("topic", colnames(tm.result$topics), sep = "")
rownames(tm.result$terms) <- paste("topic", rownames(tm.result$terms), sep = "")



dat.umap.lda <- DoUmapAndLouvain(topics.mat = posterior(out.lda)$topics, jsettings = jsettings)

# Load GLMPCA correctoin --------------------------------------------------


load(inf.glm, v=T)

tm.result.glm <- list(topics = as.matrix(glm.out$factors), terms = t(as.matrix(glm.out$loadings)))
dat.impute.glm <- t(tm.result.glm$topics %*% tm.result.glm$terms)




dat.umap.glm <- DoUmapAndLouvain(tm.result.glm$topics, jsettings = jsettings)
dat.umap.glm$louvain.glm <- dat.umap.glm$louvain; dat.umap.glm$louvain <- NULL
dat.umap.glm <- left_join(dat.umap.glm, subset(dat.umap.lda, select = c(cell, louvain)))


m.lda <- ggplot(dat.umap.lda, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = cbPalette) + ggtitle("From LDA")
m.glm <- ggplot(dat.umap.glm, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = cbPalette) + ggtitle("GLMPCA, colored by LDA louvain")
m.glm.new.louvain <- ggplot(dat.umap.glm, aes(x = umap1, y = umap2, color = louvain.glm)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
  scale_color_manual(values = cbPalette) + ggtitle("GLMPCA new louvain")

multiplot(m.lda, m.glm, m.glm.new.louvain, cols = 3)



# Plot celltypes  ---------------------------------------------------------

# jtopics <- paste("topic", c(14, 23, 4, 6, 30, 25, 11, 3), sep = "")
# names(jtopics) <- c("Eryth", "Basophil", "BcellTopic4", "Dendritic", "BcellTopic30", "Neutrophils", "InnateLymph" "Basophil2")


jtopics <- paste("topic", c(14, 3, 23, 6, 24, 11, 2, 25), sep = "")
names(jtopics) <- c("Eryth", "Basophil", "Basophil2", "Dendritic", "Bcells", "InnateLymph", "Neutrophils", "Neutrophils2")
names.final <- c("Eryth", "Basophil", "Basophil", "Dendritic", "Bcells", "InnateLymph", "Neutrophils", "Neutrophils")
jnames.hash <- hash::hash(names(jtopics), names.final)  # allows clusters to be merged afterwards by having the same final name
# names(jtopics) <- rep("Eryth", length(jtopics))

mm.celltype.lst <- FitMixtureModelLabelCells(tm.result$topics, jtopics, jthres = 0.5, show.plots = TRUE, dat.umap.long = dat.umap.lda)

# plot output
dat.celltypes <- TidyMixtureModelOutputs(mm.celltype.lst, dedup = TRUE) %>%
  rowwise() %>%
  mutate(cluster.orig = cluster,
         cluster = jnames.hash[[cluster.orig]])

# add to existing dat
dat.umap.lda <- left_join(dat.umap.lda, dat.celltypes)
dat.umap.glm <- left_join(dat.umap.glm, dat.celltypes)

m.celltype.lda <- ggplot(dat.umap.lda, aes(x = umap1, y = umap2, color = cluster)) + geom_point() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")  + 
  ggtitle(paste(jmark, jexperi))

# Plot labels in GLM corrected dlouvain  ----------------------------------

m.celltype.glm <- ggplot(dat.umap.glm, aes(x = umap1, y = umap2, color = cluster)) + geom_point() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(paste(jmark, jexperi))

m.louvain.glm <- ggplot(dat.umap.glm, aes(x = umap1, y = umap2, color = louvain.glm)) + geom_point() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(paste(jmark, jexperi))

# fill in NAs based on louvain

check.nas <- dat.umap.glm %>%
  group_by(louvain.glm, cluster) %>%
  summarise(ncell = length(cell)) %>%
  group_by(louvain.glm) %>%
  mutate(nfrac = ncell / sum(ncell))

check.nas.best <- check.nas %>%
  group_by(louvain.glm) %>%
  filter(nfrac == max(nfrac)) %>%
  ungroup()
check.nas.filt <- check.nas.best %>%
  filter(nfrac > jthres.frac.cells)

topic2clstr <- hash::hash(as.character(check.nas.filt$louvain.glm), check.nas.filt$cluster)
  
# assign NAs if nfrac > 0.9
dat.umap.glm.fillNAs <- dat.umap.glm %>% 
  rowwise() %>%
  mutate(cluster = ifelse(is.na(cluster), AssignHash(x = as.character(louvain.glm), jhash = topic2clstr, null.fill = NA), cluster))
 
m.celltype.glm.fillNAs <- ggplot(dat.umap.glm.fillNAs, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + ggtitle(paste(jmark, jexperi, "NAs imputed"))

multiplot(m.celltype.lda, m.celltype.glm, cols = 2)

print(m.celltype.glm.fillNAs)

dev.off()


