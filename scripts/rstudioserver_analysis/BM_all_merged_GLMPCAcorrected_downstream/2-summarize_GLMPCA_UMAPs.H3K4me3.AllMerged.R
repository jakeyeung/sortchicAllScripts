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

write.plots <- TRUE

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#ff9f7d", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#eb9d01", "#7fbedf", "#009E73")

jthres.frac.cells <- 0.7  # threshold for imputing NAs to celltype 
# (e.g. NA associarted to louvain cluster that are 90% of one celltype)

# topic 9 is not convergent?
jtopics <- paste("topic", c(16, 27, 25, 7, 6, 13, 9, 2, 20, 26, 14, 18), sep = "")
names(jtopics) <- c("Eryth-Sox6", "InnateLymph", "Basophils", "Eryth-Gfi1b", "pDendritic", "Bcells", "Eryth-Cdk6", "Neutrophils", "HSCs-Msi2", "HSCs-Hlf", "HSCs-Lrp5", "Dendritic")

# jtopics <- paste("topic", c(16, 27, 25, 13, 2, 26, 18), sep = "")
# names(jtopics) <- c("Eryth-Sox6", "InnateLymph", "Basophils", "Bcells", "Neutrophils", "HSCs-Hlf", "Dendritic")


names(jtopics) <- paste(names(jtopics), jtopics, sep = "_")
names.final <- names(jtopics)

# add jtopics at the end for easy reference later?
# names(jtopics) <- paste(names(jtopics), jtopics, sep = "_")
jnames.hash <- hash::hash(names(jtopics), names.final)  # allows clusters to be merged afterwards by having the same final name


jmark <- "H3K4me3"
jexperi <- "AllMerged"

mergesize <- "1000"
nbins <- "1000"
jcovar.cname <- "ncuts.var.log2.CenteredAndScaled"
jpenalty <- 1
inf.glm <- paste0("/home/jyeung/data/from_rstudioserver/scchic/rdata_robjs/GLMPCA_outputs.KeepBestPlates2/PZ_", jmark, ".", jexperi, ".KeepBestPlates2.GLMPCA_var_correction.mergebinsize_", mergesize, ".binskeep_", nbins, ".covar_", jcovar.cname, ".penalty_", jpenalty, ".winsorize_TRUE.2020-02-11.RData")
assertthat::assert_that(file.exists(inf.glm))

jinf.tss <- "/home/jyeung/hpc/databases/gene_tss/gene_tss_winsize.50000.bed"
assertthat::assert_that(file.exists(jinf.tss))

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

pdfdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/glmpca_analyses/GLMPCA_outputs.KeepBestPlates2.celltyping"
outbase <- paste0("GLMPCA_celltyping.", jmark, ".", jexperi, ".mergesize_", mergesize, ".nbins_", nbins, ".penalty_", jpenalty, ".covar_", jcovar.cname)
pdfname <- paste0(outbase, ".pdf")
outname <- paste0(outbase, ".RData")
pdfout <- file.path(pdfdir, pdfname)
outf <- file.path(pdfdir, outname)

if (write.plots){
  pdf(pdfout, useDingbats = FALSE)
}
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

# test
# jsub <- "topic9"
# mm.celltype.lst <- FitMixtureModelLabelCells(tm.result$topics, jsub, jthres = 0.5, show.plots = TRUE, dat.umap.long = dat.umap.lda, quantile.if.fail = 0.98, fail.if.nfrac = 0.25)

mm.celltype.lst <- FitMixtureModelLabelCells(tm.result$topics, jtopics, jthres = 0.5, show.plots = TRUE, dat.umap.long = dat.umap.lda, quantile.if.fail = 0.98, fail.if.nfrac = 0.25)

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
  ggtitle(paste('LDA:', jmark, jexperi))

m.louvain.glm <- ggplot(dat.umap.glm, aes(x = umap1, y = umap2, color = louvain.glm)) + geom_point() + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle(paste('GLM:', jmark, jexperi))

# fill in NAs based on louvain

check.nas <- dat.umap.glm %>%
  filter(!is.na(cluster)) %>%
  group_by(louvain.glm, cluster) %>%
  summarise(ncell = length(cell)) %>%
  group_by(louvain.glm) %>%
  mutate(nfrac = ncell / sum(ncell))

check.nas.best <- check.nas %>%
  group_by(louvain.glm) %>%
  filter(nfrac == max(nfrac)) %>%
  ungroup()

print(check.nas.best)

check.nas.filt <- check.nas.best %>%
  filter(nfrac > jthres.frac.cells)

if (nrow(check.nas.filt) == 0){
  print("No NA clusters pass threshold... consider lowering the nfrac or removing this chunk of code")
}

topic2clstr <- hash::hash(as.character(check.nas.filt$louvain.glm), check.nas.filt$cluster)
  
# assign NAs if nfrac > 0.9
dat.umap.glm.fillNAs <- dat.umap.glm %>% 
  rowwise() %>%
  mutate(cluster = ifelse(is.na(cluster), AssignHash(x = as.character(louvain.glm), jhash = topic2clstr, null.fill = NA), cluster))
 
m.celltype.glm.fillNAs <- ggplot(dat.umap.glm.fillNAs, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  ggtitle(paste('GLM:', jmark, jexperi, "NAs imputed"))


print(m.celltype.lda)
print(m.celltype.glm)
multiplot(m.celltype.lda, m.celltype.glm, cols = 2)
print(m.celltype.glm.fillNAs)


# show across conditinos 
dat.umap.glm.fillNAs <- dat.umap.glm.fillNAs %>%
  rowwise() %>%
  mutate(plate = ClipLast(cell, jsep = "-"),
         cond = GetCondFromSamp(cell, mark = jmark))
dat.umap.glm.fillNAs$cond <- factor(dat.umap.glm.fillNAs$cond, levels = c("Unenriched", "Linneg", "StemCell"))

m.celltype.glm.fillNAs.plates <- ggplot(dat.umap.glm.fillNAs, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  ggtitle(paste('GLM:', jmark, jexperi, "NAs imputed"))  + facet_wrap(~plate)

m.celltype.glm.fillNAs.conds <- ggplot(dat.umap.glm.fillNAs, aes(x = umap1, y = umap2, color = cluster)) + 
  geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  scale_color_manual(values = cbPalette, na.value = "grey85") + 
  ggtitle(paste('GLM:', jmark, jexperi, "NAs imputed"))  + facet_wrap(~cond)

print(m.celltype.glm.fillNAs.plates)
print(m.celltype.glm.fillNAs.conds)

m.celltype.glm.fillNAs.conds.nocluster <- ggplot(dat.umap.glm.fillNAs, aes(x = umap1, y = umap2, color = cond)) + 
  geom_point(alpha = 0.4) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  scale_color_manual(values = cbPalette, na.value = "grey85") +
  ggtitle(paste('GLM:', jmark, jexperi, "NAs imputed")) 
print(m.celltype.glm.fillNAs.conds.nocluster)


if (write.plots){
  dev.off()
}

# write objects to .RData
save(dat.umap.glm.fillNAs, dat.umap.glm, dat.umap.lda, mm.celltype.lst, file = outf)


