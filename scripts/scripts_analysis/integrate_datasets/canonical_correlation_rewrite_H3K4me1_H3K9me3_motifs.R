# Jake Yeung
# Date of Creation: 2019-03-18
# File: ~/projects/scchic/scripts/scripts_analysis/integrate_datasets/canonical_correlation_rewrite_H3K4me1_H3K9me3_motifs.R
# On motifs

rm(list=ls())

setwd("~/projects/scchic")

library(GGally)
library(purrr)

library(ggplot2)
library(ggrepel)
library(tidyr)
library(umap)
library(data.table)
library(dplyr)
library(hash)
library(JFuncs)
library(topicmodels)
library(scales)

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

library(cowplot)



library(GGally)

# use Seurat v3
# devtools::install_github(repo = "satijalab/seurat", ref = "release/3.0")
library(Seurat)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/IntegrateData.R")



# Get dirs ----------------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3")
mdirs <- lapply(jmarks, function(jmark){
  mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", 
                 jmark, 
                 ".filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0-",
                 "/",
                 "hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE")
  assertthat::assert_that(dir.exists(mdir))
  return(mdir)
})

jmarks.repress <- c("H3K27me3", "H3K9me3")
jbin <- "TRUE"
suffixE<- ""
suffixE1 <- ifelse(suffixE != "", paste0("--", suffixE), ".bugfix-")
suffixE2 <- ifelse(suffixE != "", paste0("_", suffixE), "")
mdirs.repress <- lapply(jmarks.repress, function(jmark){
  mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_", jbin, "-BM_", 
                 jmark, 
                 ".filt_0.99.center_TRUE", suffixE2, "-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0", suffixE1, 
                 "/",
                 "hiddenDomains_cellmin_100-cellmax_500000-binarize_", jbin, "-BM_", jmark, ".filt_0.99.center_TRUE", suffixE2)
  print(mdir)
  assertthat::assert_that(dir.exists(mdir))
  return(mdir)
})
jmarks.all <- c(jmarks, jmarks.repress)
names(jmarks.all) <- jmarks.all
mdirs <- c(mdirs, mdirs.repress)
names(mdirs) <- c(jmarks, jmarks.repress)
mara.outs <- lapply(mdirs, LoadMARA)
head(mara.outs[[1]]$act.long)


# CCA analysis ------------------------------------------------------------

jmark1 <- "H3K4me1"
jmark2 <- "H3K9me3"
# do on activation marks
X <- as.matrix(mara.outs[[jmark1]]$act.mat %>% dplyr::select(-motif)); rownames(X) <- mara.outs[[jmark1]]$act.mat$motif
# Y <- as.matrix(mara.outs$H3K4me3$act.mat %>% dplyr::select(-motif)); rownames(Y) <- mara.outs$H3K4me3$act.mat$motif
Y <- as.matrix(mara.outs[[jmark2]]$act.mat %>% dplyr::select(-motif)); rownames(Y) <- mara.outs[[jmark2]]$act.mat$motif
# Z <- cancor(X, Y, xcenter = TRUE, ycenter = TRUE)

# need to do penalized because X and Y are low rank
# source("~/projects/rmscca/scca_CVperm.R")
# source("~/projects/rmscca/scca_function.R")
# source("~/projects/rmscca/sample_sigma12_function.R")
# library(PMA)
# 
# dat <- list(X = X, Y = Y)
# # out <- PMA::MultiCCA(dat, ncomponents = 2)
# out <- PMA::CCA(x = X, z = Y, K = 2)
# 
# plot(out$u[, 1], out$u[, 2])
# plot(out$v[, 1], out$v[, 2])
# 
# X.cca <- X %*% out$u
# Y.cca <- Y %*% out$v
# 
# par(mfrow=c(1,2), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# plot(X.cca[, 1], X.cca[, 2])
# text(X.cca[, 1], X.cca[, 2], labels = mara.outs[[jmark1]]$act.mat$motif)
# plot(Y.cca[, 1], Y.cca[, 2])
# text(Y.cca[, 1], Y.cca[, 2], labels = mara.outs[[jmark2]]$act.mat$motif)

# which cells show correlation?

# plot cells with non-zero weights onto the UMAP



# plot(out$ws[[1]][, 1], out$ws[[1]][, 2])
# plot(out$ws[[2]][, 1], out$ws[[2]][, 2])


# Use Seurat implementation -----------------------------------------------

# https://satijalab.org/seurat/immune_alignment.html


# Seurat 3.0 implementation -----------------------------------------------


library(Seurat)

jmeta1 <- data.frame(cell = colnames(X), mark = "H3Km4e1", stringsAsFactors = FALSE)
jmeta3 <- data.frame(cell = colnames(Y), mark = "H3Km4e3", stringsAsFactors = FALSE)
k4me1 <- CreateSeuratObject(counts = X, meta.data = jmeta1, assay = "scChIC")
k4me3 <- CreateSeuratObject(counts = Y, meta.data = jmeta3, assay = "scChIC")

# k4me1 <- NormalizeData(k4me1, normalization.method = "RC")
k4me1 <- ScaleData(k4me1, do.scale = FALSE)
# k4me3 <- NormalizeData(k4me3, normalization.method = "RC")
k4me3 <- ScaleData(k4me3, do.scale = FALSE)

k4me1 <- FindVariableFeatures(k4me1, nfeatures = 1000)
k4me3 <- FindVariableFeatures(k4me3, nfeatures = 1000)

k4me1@meta.data$mark <- "H3K4me1"
k4me3@meta.data$mark <- "H3K4me3"

genes.use <- intersect(mara.outs$H3K4me3$act.mat$motif, mara.outs$H3K4me1$act.mat$motif)

marks.combined <- RunCCA(k4me1, k4me3, genes.use = genes.use, num.cc = 30)

p1 <- DimPlot(object = marks.combined, group.by = "mark", reduction = "cca", pt.size = 0.5, dims = 2:3)
p2 <- VlnPlot(object = marks.combined, features = "CC_1", group.by = "mark")
plot_grid(p1, p2)


# Find correlated genes ---------------------------------------------------

loadings <- marks.combined@reductions$cca@feature.loadings
embeddings <- marks.combined@reductions$cca@cell.embeddings

# reproduce cell embeddings: should be diagonalized SVD of CCA

X <- k4me1@assays$scChIC@scale.data
Y <- k4me3@assays$scChIC@scale.data

# use your own X and Y??


cca.results <- jCanonCor(X, Y, k = 20, l2.norm = TRUE)

cca.data <- rbind(cca.results$u, cca.results$v)

colnames(x = cca.data) <- paste0("CC", 1:20)
rownames(cca.data) <- c(colnames(X), colnames(Y))

# why does this matter? Maybe it does? 
cca.data.flip <- apply(cca.data, MARGIN = 2, function(x){
  if(sign(x[1]) == -1) {
    x <- x * -1
  }
  return(x)
})

# get minimum absolute number, return with actual sign
Vectorize(SelectAbsMin <- function(x1, x2){
  return(c(x1, x2)[[which.min(c(x1, x2))]])
}, vectorize.args = c("x1", "x2"), SIMPLIFY = FALSE, USE.NAMES = TRUE)


# how are these feature loadings calculated? Just projection
embeds <- cca.data.flip
# average across embeds
loads <- t(t(embeds) %*% rbind(t(X), t(Y)))
# take min function?
embeds1 <- cca.results$u
embeds2 <- cca.results$v

loads1 <- t(t(embeds1) %*% t(X))
loads2 <- t(t(embeds2) %*% t(Y))

loads.min <- matrix(mapply(SelectAbsMin, loads1, loads2), nrow = nrow(loads1), ncol = ncol(loads1), dimnames = list(rownames(loads1), colnames(loads1)))


# compare embeds and loads
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(marks.combined@reductions$cca@feature.loadings[, 1], marks.combined@reductions$cca@feature.loadings[, 2], pch = 20)
text(marks.combined@reductions$cca@feature.loadings[, 1], marks.combined@reductions$cca@feature.loadings[, 2], 
     labels = rownames(marks.combined@reductions$cca@feature.loadings))

cca1 <- 1
cca2 <- 2
plot(loads[, cca1], loads[, cca2], pch = 20)
text(loads[, cca1], loads[, cca2], 
     labels = rownames(loads))
abline(v = 0, h = 0)
# 
# plot(loads.min[, 1], loads.min[, 2], pch = 20)
# text(loads.min[, 1], loads.min[, 2], 
#      labels = rownames(loads.min))
# abline(v = 0, h = 0)


# Show zscore in both -----------------------------------------------------

act.long.merged <- rbind(mara.outs[[jmark1]]$act.long, mara.outs[[jmark2]]$act.long)
zscores.merged <- purrr::reduce(list(mara.outs[[jmark1]]$zscores, mara.outs[[jmark2]]$zscores), left_join, by = "motif")

cnames <- c("motif", paste("zscore", c(jmark1, jmark2), sep = "."))
colnames(zscores.merged) <- cnames

zscore.thres <- 0.75
zscores.merged$motif.lab <- apply(zscores.merged, 1, function(jrow){
  ifelse(max(jrow[[2]], jrow[[3]]) > zscore.thres, jrow[[1]], NA)
})


# pairs plot??
m.zscore <- ggpairs(zscores.merged, columns = paste("zscore", c(jmark1, jmark2), sep = "."), 
                    lower = list(continuous = wrap("points", alpha = 0.2))) + theme_classic()

m.zscore.both <- ggplot(zscores.merged, aes_string(x = paste0("zscore.", jmark1), y = paste0("zscore.", jmark2), label = "motif.lab")) + 
  geom_point() + geom_text_repel()

# Plot hits ---------------------------------------------------------------

nn.vec <- c(40, 40)
jmindist.vec <- c(0.2, 0.2)
jmetric <- "euclidean"
jmindist=0.1
jseed=123
custom.settings.lst <- lapply(nn.vec, function(nn) GetUmapSettings(nn, jmetric, jmindist, jseed))

meanfilt <- 10

Kstr.bin <- "15_20_25_30_35"
Kstr.nobin <- "15_20_25_30"

infs.nobin <- lapply(jmarks.all, function(jmark){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
                ".cellmin_100.cellmax_500000.binarize.", "FALSE", ".no_filt",
                "/lda_out_meanfilt.BM-", jmark, 
                ".CountThres0.K-", Kstr.nobin, ".Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})
infs.bin <- lapply(jmarks.all, function(jmark){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
                ".cellmin_100.cellmax_500000.binarize.", "TRUE", ".no_filt",
                "/lda_out_meanfilt.BM-", jmark, 
                ".CountThres0.K-", Kstr.bin, ".Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

infs <- c(infs.bin[c("H3K4me1", "H3K4me3")], infs.nobin[c("H3K27me3", "H3K9me3")])

out.lda.lst <- lapply(infs[c(jmark1, jmark2)], LoadLDA)

umap.lda.lst <- mapply(DoUmapFromLDA, out.lda.lst, custom.settings.lst, SIMPLIFY = FALSE)
names(umap.lda.lst) <- c(jmark1, jmark2)
names(mara.outs) <- jmarks.all
dat.merged.lst <- lapply(c(jmark1, jmark2), function(jmark) left_join(umap.lda.lst[[jmark]], mara.outs[[jmark]]$act.long))
names(dat.merged.lst) <- c(jmark1, jmark2)

jcolvec <- c("gray95", "gray50", "darkblue")

jmotif <- "Gata3"
jmotif <- "Cebpb"

jmotif <- "Hoxc6"

jmotif <- "Gfi1"
jmotif <- "Hmbox1"
jmotif <- "Hmbox1"

jmotif <- "Wrnip1"

jmotif <- "Pou2f2"

jmotif <- "Tal1"
jmotif <- "Cebpb"

jmotif <- "Zbtb16"

jmotif <- "Hoxc6"

m1 <- PlotMotifInUmap(jmotif, dat.merged.lst[[jmark1]], mara.outs[[jmark1]]$zscores, jmark1, jsize = 0.75, colvec = jcolvec)
m2 <- PlotMotifInUmap(jmotif, dat.merged.lst[[jmark2]], mara.outs[[jmark2]]$zscores, jmark2, jsize = 0.75, colvec = jcolvec)
multiplot(m1, m2)

# plot(loads.min[, 1], loads.min[, 2], pch = 20)
# text(loads.min[, 1], loads.min[, 2],
#      labels = rownames(loads.min))

par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0, pty = "s")
plot(loads[, 1], loads[, 2], pch = 20)
text(loads[, 1], loads[, 2],
     labels = rownames(loads))
abline(v = 0, h = 0)

# highlihgt one motif
jmot <- "Ebf1"
plot(loads[, 1], loads[, 2], pch = 20)
text(loads[, 1], loads[, 2],
     labels = sapply(rownames(loads), function(x) ifelse(x == jmot, x, "")))
abline(v = 0, h = 0)

print(m.zscore.both)

distfilt <- 0.015
# plot with radius 
loads.long <- data.frame(load1 = loads[, 1], load2 = loads[, 2], motif = rownames(loads), stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(dist = sqrt(load1 ^ 2 + load2 ^ 2)) %>%
  mutate(motif.lab = ifelse(dist > distfilt, motif, NA))
m.cca <- ggplot(loads.long, aes(x = load1, y = load2, label = motif.lab)) + geom_point() + geom_text_repel() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)
print(m.cca)

# Make plots from zscore merged and CCA -----------------------------------

jmotifs <- unlist(loads.long %>% arrange(load1) %>% filter(!is.na(motif.lab)) %>% dplyr::select(motif.lab), use.names = FALSE)

pdf(paste0("~/Dropbox/scCHiC_figs/FIG4_BM/motif_analysis/mara/", jmark1, "-", jmark2, "_CCA_analysis.", Sys.Date(), ".pdf"), useDingbats = FALSE)

print(m.zscore.both)
print(m.cca)
for (jmotif in jmotifs){
  print(jmotif)
  m1 <- PlotMotifInUmap(jmotif, dat.merged.lst[[jmark1]], mara.outs[[jmark1]]$zscores, jmark1, jsize = 0.75, colvec = jcolvec)
  m2 <- PlotMotifInUmap(jmotif, dat.merged.lst[[jmark2]], mara.outs[[jmark2]]$zscores, jmark2, jsize = 0.75, colvec = jcolvec)
  multiplot(m1, m2)
}

dev.off()

