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

FitMixtureModelLabelCells2 <-function(topics.mat, topics.keep, jthres = 0.5, show.plots = TRUE, dat.umap.long = NULL, quantile.if.fail = 0.99, fail.if.nfrac = 0.5){
  # topics.mat: rows are samples, columns are latent variables / topics
  # dat.umap.long requires columns umap1, umap2 and cell name for plotting mixturre fitting output on umap
  assertthat::assert_that(all(topics.keep %in% colnames(topics.mat)))
  
  if (is.null(names(topics.keep))){
    names(topics.keep) <- topics.keep
  }
  # label cells by topic 
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  mm.celltype.lst <- lapply(topics.keep, function(jtopic){
    print(jtopic)
    tvec.raw <- sort(topics.mat[, jtopic])
    # transform
    tvec <- log(tvec.raw / (1 - tvec.raw))
    # xline <- quantile(tvec, )
    mm <- mixtools::normalmixEM(x = tvec, lambda = c(0.9, 0.1), mu = c(-5, -1), sigma = c(2, 1), k = 2)
    (xline <- min(mm$x[which(mm$posterior[, 1] < jthres)]))
    # xline <- -2
    xcells <- names(tvec)[which(tvec > xline)]
    print(paste(length(xcells), "/", length(tvec), "assigned to", jtopic))
    if (length(xcells) / length(tvec) >= fail.if.nfrac){
      print(paste0("MM failed: too many cells assigned in mixture model, resorting to taking top: ", quantile.if.fail))
      xcells <- names(tvec)[which(tvec > quantile(tvec, quantile.if.fail))]
      print(paste("Manual threshold:", length(xcells), "/", length(tvec), "assigned to", jtopic))
    }
    if (length(xcells) == 0){
      print(paste0("No cells assigned, resorting to taking top: ", quantile.if.fail))
      xcells <- names(tvec)[which(tvec > quantile(tvec, quantile.if.fail))]
      print(paste("Manual threshold:", length(xcells), "/", length(tvec), "assigned to", jtopic))
    }
    # get topic value for assigned cells
    tvec.raw.filt <- tvec.raw[xcells]
    
    # xline <- max(mm$x[indx.btwn][which(post.filt[, 2] > 0.5)])
    if (show.plots){
      plot(density(tvec), main = paste(jtopic, jthres), xlab = "Log Odds [log(p / (1 - p))]")
      abline(v = xline, col = 'blue')
      plot.mixEM(mm, whichplots = 2, xlab2 = "Log Odds [log(p / (1 - p))]", main2 = paste(jtopic, jthres))
      abline(v = xline, col = 'blue')
    }
    
    cells.keep <- xcells
    if (show.plots){
      m.check <- PlotXYWithColor(dat.umap.long %>% mutate(is.celltype = cell %in% cells.keep), xvar = "umap1", yvar = "umap2", cname = "is.celltype", jtitle = paste(jtopic, jthres), cont.color = FALSE, col.palette = cbPalette)
      print(m.check)
    }
    
    return(list(topic = jtopic, topic.weight = tvec.raw.filt, celltype = xcells, mm = mm, threshold = xline))
  })
}


# Constants ---------------------------------------------------------------

write.plots <- TRUE

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#ff9f7d", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#eb9d01", "#7fbedf", "#009E73")

jthres.frac.cells <- 0.7  # threshold for imputing NAs to celltype 
# (e.g. NA associarted to louvain cluster that are 90% of one celltype)

# topic 9 is not convergent?
jtopics <- paste("topic", c(16, 27, 25, 7, 6, 13, 9, 2, 20, 26, 14, 18), sep = "")
names(jtopics) <- c("Eryth-Sox6", "InnateLymph", "Basophils", "Eryth-Gfi1b", "pDendritic", "Bcells", "Eryth-Cdk6", "Neutrophils", "HSCs-Msi2", "HSCs-Hlf", "HSCs-Lrp5", "Dendritic")

# jtopics <- paste("topic", c(16, 27, 25, 7, 6, 13, 2, 20, 26, 14, 18), sep = "")
# names(jtopics) <- c("Eryth-Sox6", "InnateLymph", "Basophils", "Eryth-Gfi1b", "pDendritic", "Bcells", "Neutrophils", "HSCs-Msi2", "HSCs-Hlf", "HSCs-Lrp5", "Dendritic")

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
pdfname <- paste0("GLMPCA_celltyping.", jmark, ".", jexperi, ".mergesize_", mergesize, ".nbins_", nbins, ".penalty_", jpenalty, ".covar_", jcovar.cname, ".pdf")
pdfout <- file.path(pdfdir, pdfname)

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

if (write.plots){
  dev.off()
}


