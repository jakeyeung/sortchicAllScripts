# Jake Yeung
# Date of Creation: 2019-02-05
# File: ~/projects/scchic/scripts/scripts_analysis/debug_marker_genes/color_H3K27me3_by_metacell.R
# Bin analysis, color by metacell to figure out what's going on

rm(list=ls())

library(dplyr)
library(ggplot2)
library(JFuncs)
library(topicmodels)
library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(hash)
library(umap)
library(metacell)
library(scales)

source("scripts/Rfunctions/PlotFunctions.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")

jmark <- "H3K27me3"
jbin <- TRUE
if (jbin){
  Kstr <- "5_10_15_20_25"
} else {
  Kstr <- "5_15_25"
}
inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.", 
              jbin, "/lda_out_meanfilt.BM-", 
              jmark, ".CountThres0.K-", Kstr, ".Robj")
inf.nobin <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.", 
              FALSE, "/lda_out_meanfilt.BM-", 
              jmark, ".CountThres0.K-", "5_15_25", ".Robj")
inf.mc <- file.path(paste0("~/Dropbox/scCHiC_figs/FIG4_BM/", jmark, ".datadir_mc_f.Rda"))

assertthat::assert_that(file.exists(inf))
assertthat::assert_that(file.exists(inf.nobin))
assertthat::assert_that(file.exists(inf.mc))


# Load  -------------------------------------------------------------------

load(inf.nobin, v=T)
count.mat.nobin <- count.mat
load(inf, v=T)
load(inf.mc, v=T)

out.lda <- ChooseBestLDA(out.lda)
kchoose <- out.lda@k

tm.result <- posterior(out.lda)
top.thres <- 0.96 
topic.regions <- lapply(seq(kchoose), function(clst){
  return(SelectTopRegions(tm.result$terms[clst, ], colnames(tm.result$terms), method = "thres", method.val = top.thres))
})
# cluster the terms, but first get interesting terms
top.regions <- unique(unlist(topic.regions))  # topic.regions defined by threshold, 98 or 96th percentile of top weights in each column of the betas matrix

terms.mat <- t(tm.result$terms)[top.regions, ]

regions <- data.frame(seqnames = sapply(colnames(tm.result$terms), GetChromo),
                      start = sapply(colnames(tm.result$terms), GetStart),
                      end = sapply(colnames(tm.result$terms), GetEnd), 
                      stringsAsFactors = FALSE)
rownames(regions) <- colnames(tm.result$terms)
regions <- subset(regions, !seqnames %in% c("chr20", "chr21"))

regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
                                                annoDb='org.Mm.eg.db'))
regions.annotated$region_coord <- names(regions.range)


# Set up metacell ---------------------------------------------------------

mc.out <- GetMCObjects(inf = inf.mc)
mc_index <- mc.out$mc_index; mc_colors <- mc.out$mc_colors

cellnames <- unname(out.lda@documents)

cells.clstr.hash <- hash(names(mc_index), mc_index)
cells.color.hash <- hash(as.character(seq(mc_colors)), mc_colors)

# filter for common cells between the two analyses
topics.mat <- tm.result$topics
cells.lda <- unname(rownames(topics.mat))
cells.metacell <- names(mc_index)

cells.common <- intersect(cells.lda, cells.metacell)
print(paste("N cells:", length(cells.common)))

topics.mat <- topics.mat[cells.common, ]
cells.indx <- sapply(rownames(topics.mat), function(x) cells.clstr.hash[[x]])
cells.rgb <- sapply(cells.indx, function(x) cells.color.hash[[as.character(x)]])



nn=40
jmetric='euclidean' 
jmindist=0.05
custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist)

dat.umap <- umap(topics.mat, config = custom.settings)
rownames(dat.umap$layout) <- rownames(topics.mat)
jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)

# compare metacell and lda
plot(dat.umap$layout[, 1], dat.umap$layout[, 2], 
     pch = 20, 
     main = paste0("MC colors: ", jmark), 
     col = cells.rgb, 
     asp = 1, 
     xlab = "UMAP Dim 1", ylab = "UMAP Dim 2")



# Too few cells? ----------------------------------------------------------




# Show Hox counts ---------------------------------------------------------


# [1] 17
# [1] "chr15:103040000-103140000"
# [1] "chr15:103000000-103100000"
# [1] "chr15:102980000-103080000"
# [1] "chr15:102960000-103060000"
# [1] 17
# [1] "chr15:102840000-102940000"
# mat.norm <- as.matrix(count.mat.nobin)
mat.norm <- as.matrix(count.mat)
mat.norm <- log10(sweep(mat.norm, MARGIN = 2, STATS = colSums(mat.norm), FUN = "/") + 10^-6)
mat.imputed <- t(tm.result$topics %*% tm.result$terms)

# mat.norm <- log2(sweep(mat.norm, MARGIN = 2, STATS = colSums(mat.norm), FUN = "/") * 10^6 + 1)
# plot(density(mat.norm))
# mat.norm <- log2(mat.norm + 1)


jpeaks <- c("chr15:103040000-103140000", "chr15:103000000-103100000", "chr15:102980000-103080000", 
            "chr15:102960000-103060000", "chr15:102840000-102940000")

pdf("~/Dropbox/scCHiC_figs/FIG4_BM/marker_genes/H3K27me3_topic23_top_hits.pdf")
for (jpeak in jpeaks){
  jcounts.norm <- mat.norm[jpeak, ]
  # imputed counts
  jcounts.imputed <- log10(mat.imputed[jpeak, ])
  
  dat <- data.frame(umap1 = dat.umap$layout[, 1], 
                    umap2 = dat.umap$layout[, 2], 
                    counts.norm = jcounts.norm,
                    counts.imputed = jcounts.imputed)
  print(head(dat))
  jsize <- 0.5
  m.count <- ggplot(dat, aes(x = umap1, y = umap2, col = counts.norm)) + 
    geom_point(size = jsize) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size=7), legend.position = "bottom") + 
    scale_color_gradient2(low = muted("blue"), mid = "white", high = muted("red"), 
                          midpoint = mean(dat$counts.norm)) + 
    ggtitle(jpeak)
  m.imputed <- ggplot(dat, aes(x = umap1, y = umap2, col = counts.imputed)) + 
    geom_point(size = jsize) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size=7), legend.position = "bottom") + 
    scale_color_gradient2(low = muted("blue"), mid = "white", high = muted("red"), 
                          midpoint = mean(dat$counts.imputed)) + 
    ggtitle(jpeak)
  multiplot(m.count, m.imputed, cols = 2)
}

for (i in seq(10)){
  jpeak <- names(sort(tm.result$terms[23, ], decreasing = TRUE))[[i]]

  jcounts.norm <- mat.norm[jpeak, ]
  
  # imputed counts
  mat.imputed <- t(tm.result$topics %*% tm.result$terms)
  # jcounts.imputed <- log2(mat.imputed[jpeak, ] * 10^6 + 1)
  jcounts.imputed <- log10(mat.imputed[jpeak, ])
  
  dat <- data.frame(umap1 = dat.umap$layout[, 1], 
                    umap2 = dat.umap$layout[, 2], 
                    counts.norm = jcounts.norm,
                    counts.imputed = jcounts.imputed)
  print(head(dat))
  jsize <- 0.5
  m.count <- ggplot(dat, aes(x = umap1, y = umap2, col = counts.norm)) + 
    geom_point(size = jsize) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size=7), legend.position = "bottom") + 
    scale_color_gradient2(low = muted("blue"), mid = "white", high = muted("red"), 
                          midpoint = mean(dat$counts.norm)) + 
    ggtitle(jpeak)
  # print(m.count)
  
  m.imputed <- ggplot(dat, aes(x = umap1, y = umap2, col = counts.imputed)) + 
    geom_point(size = jsize) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size=7), legend.position = "bottom") + 
    scale_color_gradient2(low = muted("blue"), mid = "white", high = muted("red"), 
                          midpoint = mean(dat$counts.imputed)) + 
    ggtitle(jpeak)
  # print(m.imputed)
  multiplot(m.count, m.imputed, cols = 2)
  
}
  
dev.off()
  
  
