# Jake Yeung
# Date of Creation: 2019-03-11
# File: ~/projects/scchic/scripts/scripts_analysis/integrate_datasets/multi_canonical_correlation_rewrite_redo.R
# Redo with own implementation 

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

library(irlba)

library(GGally)

# use Seurat v3
# devtools::install_github(repo = "satijalab/seurat", ref = "release/3.0")
library(Seurat)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")




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
mdirs.repress <- lapply(jmarks.repress, function(jmark){
  mdir <- paste0("/Users/yeung/data/scchic/from_cluster/mara_analysis/hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", 
                 jmark, 
                 ".filt_0.99.center_TRUE-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix-",
                 "/",
                 "hiddenDomains_cellmin_100-cellmax_500000-binarize_TRUE-BM_", jmark, ".filt_0.99.center_TRUE")
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
jmark2 <- "H3K4me3"
jmark3 <- "H3K27me3"
jmark4 <- "H3K9me3"

# do on activation marks
X <- as.matrix(mara.outs[[jmark1]]$act.mat %>% dplyr::select(-motif)); rownames(X) <- mara.outs[[jmark1]]$act.mat$motif
Y <- as.matrix(mara.outs[[jmark2]]$act.mat %>% dplyr::select(-motif)); rownames(Y) <- mara.outs[[jmark2]]$act.mat$motif
Z <- as.matrix(mara.outs[[jmark3]]$act.mat %>% dplyr::select(-motif)); rownames(Z) <- mara.outs[[jmark3]]$act.mat$motif
A <- as.matrix(mara.outs[[jmark4]]$act.mat %>% dplyr::select(-motif)); rownames(A) <- mara.outs[[jmark4]]$act.mat$motif

dat <- list(X = X, Y = Y, Z = Z, A = A)

out <- jMultiCCA(dat, num.ccs = 10, niter = 25)

ws.norm <- apply(out$ws[[1]], 2, L2Norm)



# analyze output???
embeds.lst <- out$ws
loadings.lst <- mapply(function(W, X) return(X %*% W), out$ws, dat, SIMPLIFY = FALSE)

# take average?? or min??
loadings.min <- matrix(mapply(SelectAbsMin2, loadings.lst[[1]], loadings.lst[[2]], loadings.lst[[3]], loadings.lst[[4]], SIMPLIFY = FALSE, USE.NAMES = TRUE), 
                       nrow = nrow(loadings.lst[[1]]), ncol = ncol(loadings.lst[[1]]), 
                       dimnames = list(rownames(loadings.lst[[1]]), colnames(loadings.lst[[2]])))

loadings.3d <- array(data = unlist(loadings.lst), c(677, 10, 4), dimnames = list((rownames(loadings.lst[[1]]))))
loadings.mean <- apply( loadings.3d , 1:2 , mean ); rownames(loadings.mean) <- rownames(loadings.3d)
loadings.median <- apply( loadings.3d , 1:2 , median ); rownames(loadings.median) <- rownames(loadings.3d)
# loadings.mean <- matrix(mapply(mean, loadings.lst[[1]], loadings.lst[[2]], loadings.lst[[3]], loadings.lst[[4]], SIMPLIFY = FALSE, USE.NAMES = TRUE), 
#                        nrow = nrow(loadings.lst[[1]]), ncol = ncol(loadings.lst[[1]]), 
#                        dimnames = list(rownames(loadings.lst[[1]]), colnames(loadings.lst[[2]])))
# 
# loadings.median <- matrix(mapply(median, loadings.lst[[1]], loadings.lst[[2]], loadings.lst[[3]], loadings.lst[[4]], SIMPLIFY = FALSE, USE.NAMES = TRUE), 
#                        nrow = nrow(loadings.lst[[1]]), ncol = ncol(loadings.lst[[1]]), 
#                        dimnames = list(rownames(loadings.lst[[1]]), colnames(loadings.lst[[2]])))

cc1 <- 1
cc2 <- 2

cc1 <- 1
cc2 <- 2
plot(loadings.min[, cc1], loadings.min[, cc2], pch = 20, main = paste(cc1, "vs", cc2))
text(loadings.min[, cc1], loadings.min[, cc2], labels = rownames(loadings.min))

cc1 <- 5
cc2 <- 6
plot(loadings.mean[, cc1], loadings.mean[, cc2], pch = 20, main = paste(cc1, "vs", cc2))
text(loadings.mean[, cc1], loadings.mean[, cc2], labels = rownames(loadings.mean))

plot(loadings.median[, cc1], loadings.median[, cc2], pch = 20, main = paste(cc1, "vs", cc2))
text(loadings.median[, cc1], loadings.median[, cc2], labels = rownames(loadings.median))


# Plot examples  ----------------------------------------------------------

nn.vec <- c(40, 35, 40, 40)
jmindist.vec <- c(0.2, 0.1, 0.2, 0.1)
jmetric <- "euclidean"
jseed=123
# custom.settings.lst <- lapply(nn.vec, function(nn) GetUmapSettings(nn, jmetric, jmindist, jseed))
custom.settings.lst <- mapply(function(nn, jmindist) GetUmapSettings(nn, jmetric, jmindist, jseed), nn.vec, jmindist.vec, SIMPLIFY = FALSE)

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
# out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks.all, infs, SIMPLIFY = FALSE)

out.lda.lst <- lapply(infs, LoadLDA)

# out.lda.lst <- out.objs

umap.lda.lst <- mapply(DoUmapFromLDA, out.lda.lst, custom.settings.lst, SIMPLIFY = FALSE)

names(umap.lda.lst) <- jmarks.all
names(mara.outs) <- jmarks.all

dat.merged.lst <- lapply(jmarks.all, function(jmark) left_join(umap.lda.lst[[jmark]], mara.outs[[jmark]]$act.long))
names(dat.merged.lst) <- jmarks.all

jmotif <- "Hoxc6"

jmotif <- "Cebpb"

jmotif <- "Gata3"
jmotif <- "Gfi1"
jmotif <- "Hoxa2"
jmotif <- "Bcl3"

jmotif <- "Tal1"
jmotif <- "Zbtb16"
jmotif <- "Hoxa2"
jmotif <- "Hoxa1"
jmotif <- "Zfp110"
jmotif <- "Gfi1"
jmotif <- "Foxc1"


jmotif <- "Hoxb6"

jmotif <- "Ebf1"
jmotif <- "Ezh2"

jmotif <- "Ezh2"

jmotif <- "Ebf1"

jmotif <- "Nfatc1"
jmotif <- "Hoxa1"
jmotif <- "Zbtb16"

jmotif <- "Stat6"
m1 <- PlotMotifInUmap(jmotif, dat.merged.lst[[jmark1]], mara.outs[[jmark1]]$zscores, jmark1, jsize = 0.75)
m2 <- PlotMotifInUmap(jmotif, dat.merged.lst[[jmark2]], mara.outs[[jmark2]]$zscores, jmark2, jsize = 0.75)
m3 <- PlotMotifInUmap(jmotif, dat.merged.lst[[jmark3]], mara.outs[[jmark3]]$zscores, jmark3, jsize = 0.75)
m4 <- PlotMotifInUmap(jmotif, dat.merged.lst[[jmark4]], mara.outs[[jmark4]]$zscores, jmark4, jsize = 0.75)
multiplot(m1, m2, m3, m4, cols = 2)

cc1 <- 1
cc2 <- 2
plot(loadings.min[, cc1], loadings.min[, cc2], pch = 20, main = paste(cc1, "vs", cc2))
text(loadings.min[, cc1], loadings.min[, cc2], labels = rownames(loadings.min))

plot(loadings.mean[, cc1], loadings.min[, cc2], pch = 20, main = paste(cc1, "vs", cc2))
text(loadings.mean[, cc1], loadings.min[, cc2], labels = rownames(loadings.mean))

plot(loadings.median[, cc1], loadings.median[, cc2], pch = 20, main = paste(cc1, "vs", cc2))
text(loadings.median[, cc1], loadings.median[, cc2], labels = rownames(loadings.median))


# loads.umap <- umap(loads.min[, 1:10])
# plot(loads.umap$layout[, 1], loads.umap$layout[, 2], pch = 20)
# text(loads.umap$layout[, 1], loads.umap$layout[, 2], labels = rownames(loads.umap$layout))

# Now let’s do multiCCA ---------------------------------------------------



