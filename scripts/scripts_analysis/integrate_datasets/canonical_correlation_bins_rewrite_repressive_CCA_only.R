# Jake Yeung
# Date of Creation: 2019-03-18
# File: ~/projects/scchic/scripts/scripts_analysis/integrate_datasets/canonical_correlation_bins_rewrite_repressive_CCA_only.R
# CCA only 

# https://satijalab.org/seurat/immune_alignment.html
# https://satijalab.org/seurat/Seurat_AlignmentTutorial.html

rm(list=ls())

setwd("~/projects/scchic")

jstart <- Sys.time()


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

library(ggrastr)

# use Seurat v3
# devtools::install_github(repo = "satijalab/seurat", ref = "release/3.0")
library(Seurat)

source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/PlotFunctions.R")

source("scripts/Rfunctions/IntegrateData.R")


# Load active marks LDA output --------------------------------------------


# marks.keep <- c("H3K4me1", "H3K4me3")
marks.keep <- c("H3K27me3", "H3K9me3")

jsize <- 0.5
# jcolvec <- c("blue", "yellow", "red")
# jcolvec <- c("blue", "gray80", "red")
jcolvec <- c("gray70", "gray50", "darkblue")


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
names(jmarks) <- jmarks

meanfilt <- 10

Kstr.bin <- "15_20_25_30_35"
Kstr.nobin <- "15_20_25_30"

infs.nobin <- lapply(jmarks, function(jmark){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
                ".cellmin_100.cellmax_500000.binarize.", "FALSE", ".no_filt",
                "/lda_out_meanfilt.BM-", jmark, 
                ".CountThres0.K-", Kstr.nobin, ".Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})
infs.bin <- lapply(jmarks, function(jmark){
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_", meanfilt, 
                ".cellmin_100.cellmax_500000.binarize.", "TRUE", ".no_filt",
                "/lda_out_meanfilt.BM-", jmark, 
                ".CountThres0.K-", Kstr.bin, ".Robj")
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

infs <- c(infs.bin[c("H3K4me1", "H3K4me3")], infs.nobin[c("H3K27me3", "H3K9me3")])
out.objs <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks, infs, SIMPLIFY = FALSE)
out.objs.nobin <- mapply(function(jmark, inf) LoadLDABins(jmark, inf=inf), jmarks, infs.nobin, SIMPLIFY = FALSE)
# out.objs <- mapply(function(jmark, inf) print(paste(jmark, inf)), jmarks, infs)
names(out.objs) <- jmarks
names(out.objs.nobin) <- jmarks

tm.result.lst <- lapply(out.objs, function(x) posterior(x[['out.lda']]))

mat.impute.lst <- lapply(tm.result.lst, function(tm.result) t(tm.result$topic %*% tm.result$term))

# use nobin for mat
count.mat.lst <- lapply(out.objs.nobin, function(x) sweep(as.matrix(x$count.mat), 2, Matrix::colSums(x$count.mat), "/"))

print(lapply(mat.impute.lst, dim))

# Wrangle data ------------------------------------------------------------

top.thres <- 0.999
top.regions.lst <- lapply(jmarks, function(jmark){
  topic.regions <- out.objs[[jmark]]$topic.regions  # for each cluster
  # topic.regions <- lapply(out.objs[[jmark]]$out.lda, function(clst){
  #   return(SelectTopRegions(tm.result$terms[clst, ], colnames(tm.result$terms), method = "thres", method.val = top.thres))
  # })
  top.regions <- unique(unlist(topic.regions))
  
  top.regions <- unique(unlist(topic.regions))
  return(top.regions) 
})
top.regions.merge <- unique(unlist(top.regions.lst))

mat.merged.lst <- lapply(mat.impute.lst, function(mat.impute){
  row.i <- which(top.regions.merge %in% rownames(mat.impute))
  return(as.data.frame(mat.impute[row.i, ]))
})

rnames.lst <- lapply(mat.merged.lst, function(x) rownames(x))
top.regions.intersect <- purrr::reduce(.x = rnames.lst, .f = intersect)

print(length(top.regions.intersect))

mat.merged.lst <- lapply(mat.merged.lst, function(mat.merged){
  return(mat.merged[top.regions.intersect, ])
})


# Do penalized PMD? -------------------------------------------------------


library(Seurat)


X <- as.matrix(mat.merged.lst[[marks.keep[[1]]]])
Y <- as.matrix(mat.merged.lst[[marks.keep[[2]]]])

X <- scale(X, center = TRUE, scale = FALSE)
Y <- scale(Y, center = TRUE, scale = FALSE)

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

cca1 <- 1
cca2 <- 2
plot(loads[, cca1], loads[, cca2], pch = 20)
# text(loads[, cca1], loads[, cca2], 
#      labels = rownames(loads))
abline(v = 0, h = 0)


distfilt <- 60



loads.long <- data.frame(loads, bin = rownames(loads), stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(dist = sqrt(CC1 ^ 2 + CC2 ^ 2)) %>%
  mutate(bin.lab = ifelse(dist > distfilt, bin, NA))
# add gene info
jmark <- marks.keep[[1]]
loads.long <- left_join(loads.long, out.objs[[jmark]]$regions.annot %>% dplyr::select(c(SYMBOL, region_coord)), by = c("bin" = "region_coord"))
loads.long$bingene <- paste(loads.long$bin, loads.long$SYMBOL, sep = ";")

m.cca <- ggplot(loads.long, aes(x = CC1, y = CC2, label = bin.lab)) + geom_point() + geom_text_repel() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)
print(m.cca)


# Plot some top hits ------------------------------------------------------

# plto constants

jscale.fac <- 10^6
jpseudo <- 1

# can save time by precalculating the UMAP and feeding it into the plot functions, also can customize for each UMAP
# custom settings for each UMAP
jmetric='euclidean'
jmindist=0.2
jseed=123
nn.vec <- c(40, 35, 40, 40)
jmindist.vec <- c(0.2, 0.1, 0.2, 0.1)
custom.settings.lst <- mapply(function(nn, jmindist) GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist, seed = jseed), 
                              nn.vec, jmindist.vec, SIMPLIFY = FALSE)
# plot umaps to check
topics.mat.lst <- lapply(out.objs, function(x) x$tm.result$topics)
dat.umap.lst <- mapply(function(custom.settings, topics.mat){
  dat.umap <- umap(topics.mat, config = custom.settings) 
  return(dat.umap)
}, custom.settings.lst, topics.mat.lst, SIMPLIFY = FALSE)
names(dat.umap.lst) <- jmarks

# jmark <- "H3K4me1"
jpeak <- "chr7:73520000-73620000"  # top loads1 peak

# jpeak <- (loads.long %>% arrange(desc(load2)))$bin[[1]]  # highest load peak
# jpeak <- (loads.long %>% arrange(load2))$bin[[1]]  # lowest



jpeak <- (loads.long %>% arrange(CC3))$bin[[1]]
(jgene <- subset(out.objs[[jmark]]$regions.annot, region_coord == jpeak)$SYMBOL[[1]])

loads.long %>% arrange(CC3)
loads.long %>% arrange(desc(CC3))
loads.long %>% arrange(CC1)
loads.long %>% arrange(desc(CC1))
loads.long %>% arrange(desc(CC2))

jpeak <- "chr6:47560000-47660000"
(jgene <- subset(out.objs[[jmark]]$regions.annot, region_coord == jpeak)$SYMBOL[[1]])
system.time(
  PlotUmapAllMarks(jmarks[marks.keep], tm.result.lst[marks.keep], jpeak, juse.count.mat = count.mat.lst[marks.keep], dat.umap.lst[marks.keep], jgene, jsize, jcolvec, .log = TRUE, scale.fac = jscale.fac, pseudocount = jpseudo)
)

# plot CCA outputs
# plot top CC1s
jtopn <- 30
# loads.dat %>% arrange(desc(CC1))
# take top 30 unique genes
bins.to.plot1 <- (loads.long %>% arrange(desc(CC1)) %>% group_by(SYMBOL) %>% filter(CC1 == max(CC1)))$bingene[1:jtopn]
bins.to.plot1b <- (loads.long %>% arrange(CC1) %>% group_by(SYMBOL) %>% filter(CC1 == max(CC1)))$bingene[1:jtopn]
bins.to.plot2 <- (loads.long %>% arrange(desc(CC2)) %>% group_by(SYMBOL) %>% filter(CC2 == max(CC2)))$bingene[1:jtopn]
bins.to.plot2b <- (loads.long %>% arrange(CC2) %>% group_by(SYMBOL) %>% filter(CC2 == min(CC2)))$bingene[1:jtopn]
bins.to.plot3 <- (loads.long %>% arrange(CC3) %>% group_by(SYMBOL) %>% filter(CC3 == min(CC3)))$bingene[1:jtopn]
bins.to.plot3b <- (loads.long %>% arrange(desc(CC3)) %>% group_by(SYMBOL) %>% filter(CC3 == min(CC3)))$bingene[1:jtopn]

pdf(paste0("~/Dropbox/scCHiC_figs/FIG4_BM/analyses/2019-03-17_cca_output_", paste(marks.keep, collapse = "-"), "_CCA_only.pdf"), useDingbats = FALSE)
jsize <- 0.1
m.cca <- ggplot(loads.long %>% mutate(bin.lab = ifelse(bingene %in% bins.to.plot1, bingene, NA)) %>% filter(!is.na(bingene)), aes(x = CC1, y = CC2, label = bin.lab)) + 
  # geom_point(alpha = 0.2) + 
  geom_point_rast(alpha = 0.2, color = "blue", size = jsize) +
  geom_text_repel() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)
print(m.cca)


bins.to.plot.lst <- list(bins.to.plot1, bins.to.plot1b, bins.to.plot2, bins.to.plot2b, bins.to.plot3, bins.to.plot3b)

lapply(bins.to.plot.lst, function(bins.to.plot){
  m.cca <- ggplot(loads.long %>% mutate(bin.lab = ifelse(bingene %in% bins.to.plot, bingene, NA)), aes(x = CC1, y = CC2, label = bin.lab)) + 
    # geom_point(alpha = 0.2) + 
    geom_point_rast(alpha = 0.2, color = "pink", size = jsize) + 
    geom_text_repel() + geom_vline(xintercept = 0) + geom_hline(yintercept = 0)
  print(m.cca)
})
dev.off()
print(Sys.time() - jstart)
