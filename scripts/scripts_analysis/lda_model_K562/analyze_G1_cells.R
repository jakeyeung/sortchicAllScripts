# Jake Yeung
# Date of Creation: 2019-02-17
# File: ~/projects/scchic/scripts/scripts_analysis/lda_model_K562/analyze_G1_cells.R
# Analyze G1 cells 


rm(list=ls())

library(dplyr)
library(ggplot2)
library(topicmodels)
library(cisTopic)
library(GenomicRanges)
library(ChIPseeker)
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(hash)
library(umap)
library(tidytext)

library(igraph)  # louvain

library(Gviz)

library(biomaRt)

source("scripts/Rfunctions/ParseStrings.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/PlotFunctions.R")

# Load data ---------------------------------------------------------------

# jmark <- "H3K27me3"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")

for (jmark in jmarks){
  
  # jmark <- "H3K9me3"
  # jmark <- "H3K4me1"
  # jmark <- "H3K4me3"
  # jbin <- "TRUE"; kstr <- "15_20_25_30_35"
  jbin <- "FALSE"; kstr <- "15_20_25_30"
  dirmain <- paste0("/Users/yeung/data/scchic/from_cluster/LDA_out_K562/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_10.cellmin_100.cellmax_500000.binarize.", jbin, ".no_filt")
  inf <- file.path(dirmain, paste0("lda_out_meanfilt.K562-G1-", jmark, ".CountThres0.K-", kstr, ".Robj"))
  
  assertthat::assert_that(file.exists(inf))
  
  load(inf, v=T)
  
  out.lda <- ChooseBestLDA(out.lda)
  kchoose <- out.lda@k
  
  tm.result <- posterior(out.lda)
  
  top.thres <- 0.995
  topic.regions <- lapply(seq(kchoose), function(clst){
    return(SelectTopRegions(tm.result$terms[clst, ], colnames(tm.result$terms), method = "thres", method.val = top.thres))
  })
  top.regions <- unique(unlist(topic.regions))
  
  
  # Process regions ---------------------------------------------------------
  
  # chainpath <- "~/data/scchic/databases/hg38ToHg19.over.chain"
  # assertthat::assert_that(file.exists(chainpath))
  # ch <- rtracklayer::import.chain(chainpath)
  # seqlevelsStyle(regions.range) = "UCSC"
  
  regions <- data.frame(seqnames = sapply(colnames(tm.result$terms), GetChromo),
                        start = sapply(colnames(tm.result$terms), GetStart),
                        end = sapply(colnames(tm.result$terms), GetEnd),
                        terms = colnames(tm.result$terms), 
                        stringsAsFactors = FALSE) %>%
    filter(seqnames != "chr23" &  start != "")  # negative values are returned as "", remove them
  
  
  rownames(regions) <- regions$terms
  
  regions.range <- makeGRangesFromDataFrame(as.data.frame(regions)) 
  
  regions.annotated <- as.data.frame(annotatePeak(regions.range,
                                                  TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
                                                  annoDb='org.Hs.eg.db'))
  regions.annotated <- regions.annotated %>%
    rowwise() %>%
    mutate(region_coord = paste(paste(seqnames, start, sep = ":"), end, sep = "-"))
  # regions.annotated$hg38peak <- names(regions.range) 
  
  # regions.range.19 = unlist(rtracklayer::liftOver(regions.range, ch))
  # regions.range.19$hg38peak <- names(regions.range.19)
  # 
  # regions.annotated.19 <- as.data.frame(annotatePeak(unname(regions.range.19),
  # TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene,
  # annoDb='org.Hs.eg.db'))
  
  
  # CSet up plots -----------------------------------------------------------
  
  
  nn <- 40
  nn.terms <- 15
  jmetric <- 'euclidean'
  jmindist <- 0.4
  custom.settings <- GetUmapSettings(nn, jmetric, jmindist)
  custom.settings.terms <- GetUmapSettings(nn.terms, jmetric, jmindist)
  
  
  dat.umap <- umap(tm.result$topics, config = custom.settings)
  rownames(dat.umap$layout) <- rownames(tm.result$topics)
  
  jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)
  
  
  
  # Do plots ----------------------------------------------------------------
  
  
  jcounts <- Matrix::colSums(count.mat)
  jcol.counts <- ColorsByCounts(jcounts, nbreaks = 100)
  
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = jmain, pty = 's', col = jcol.counts, asp = 0.5)
  
  
  dat.pca <- prcomp(tm.result$topics, center = TRUE, scale. = TRUE)
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  plot(dat.pca$x[, 1], dat.pca$x[, 2], pch = 20, col = jcol.counts)
  
  jcol.rgbs <- lapply(seq(kchoose), ColorsByGamma, tm.result$topics)
  nb.col <- 5
  nb.row <- ceiling(kchoose / nb.col)
  par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
  mapply(function(jcol.rgb, jtopic){
    plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.5)
  }, jcol.rgbs, seq(kchoose))
  
  # Do plots with alpha no color, ggplot2 style 
  
  dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], stringsAsFactors = FALSE)
  
  # m1 <- ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(alpha = 0.5, size = 4)
  
  # Check top hits ----------------------------------------------------------
  
  top.peaks <- tidytext::tidy(out.lda, matrix = "beta") %>%
    group_by(topic) %>%
    arrange(desc(beta)) %>%
    mutate(rnk = seq(length(beta)))
  
  jtopic <- 28
  topn <- 50
  # find top hits
  jsub <- subset(top.peaks, topic == jtopic)
  
  jpeaks <- jsub$term[1:topn]
  
  print(jpeaks)
  
  print(head(subset(regions.annotated, region_coord %in% jpeaks)))
  
  
  
  
  # Do louvain on this ------------------------------------------------------
  
  dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout), stringsAsFactors = FALSE)
  
  cell.indx <- hash(rownames(dat.umap$knn$indexes), dat.umap$knn$indexes[, 1])
  cell.indx.rev <- hash(dat.umap$knn$indexes[, 1], rownames(dat.umap$knn$indexes))
  nr <- nrow(dat.umap$knn$indexes)
  nc <- ncol(dat.umap$knn$indexes)
  # nc <- 4
  # edgelist <- matrix(NA, nrow = nr * nc, ncol = 3)
  edgelist <- matrix(NA, nrow = nr * nc, ncol = 2)
  colnames(edgelist) <- c("from", "to")
  for (vertex.i in seq(nr)){
    istart <- nc*(vertex.i - 1)+1
    iend <- nc*vertex.i
    edgelist[istart : iend, 1] <- cell.indx.rev[[as.character(vertex.i)]]
    edgelist[istart : iend, 2] <- sapply(dat.umap$knn$indexes[vertex.i, 1:nc], function(x) cell.indx.rev[[as.character(x)]])
    # edgelist[istart : iend, 3] <- 1 / (dat.umap$knn$distances[vertex.i, 1:nc] + 0.1)
  }
  g <- graph_from_data_frame(edgelist, directed=FALSE)
  g.out <- cluster_louvain(g, weights = NULL)
  V(g)$color <- g.out$membership
  clstr <- hash(g.out$names, g.out$membership)
  
  dat.umap.long$louvain <- sapply(dat.umap.long$cell, function(x) clstr[[x]])
  
  # plot umap with louvain
  m.louvain <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = as.character(louvain))) + geom_point() +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_brewer(palette = "Spectral")
  print(m.louvain)
  
  
  
  # What are the fold changes??? --------------------------------------------
  
  mat.norm <- t(tm.result$topics %*% tm.result$terms)  # this should give normalized signal, without the poisson noise?
  
  jpeak <- top.peaks$term[[1]]
  
  jchromo <- strsplit(jpeak, ":")[[1]][[1]]
  
  # get data for chromo 13?
  jpeaks <- grep(jchromo, rownames(mat.norm), value = TRUE)
  x <- as.data.frame(mat.norm[jpeaks, ])
  x.long <- data.frame(exprs = unlist(x), cell = rep(colnames(x), each = nrow(x)),
                       coord = rep(rownames(x), ncol(x)), stringsAsFactors = FALSE)
  x.long$louvain.orig <- sapply(x.long$cell, function(x) clstr[[x]])
  # x.long$louvain <- sapply(as.character(x.long$louvain.orig), function(x) remap.clstr[[x]])
  x.long$louvain <- x.long$louvain.orig
  x.long$exprs <- x.long$exprs * 10^6
  
  # remove negative starts
  x.long <- subset(x.long, !grepl(":-", coord))
  
  
  # Plot a top hit ----------------------------------------------------------
  
  mart.obj <- useMart(biomart = 'ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl', host="www.ensembl.org")
  
  jgene <- ""
  # jpeak <- "chr13:91300000-91400000"
  # pick top peka in chr13?
  # jpeak <- subset(top.peaks, grepl("chr13", term))$term[[1]]
  # jpeak <- "chr1:143260000-143360000"
  # jgene <- "Sox6"
  # jpeak <- subset(top.peaks.annotated, topic == jtopic & SYMBOL == jgene)$term[[1]]
  
  print(PlotImputedPeaks(tm.result, jpeak, jmark, show.plot = FALSE,
                         return.plot.only = TRUE, usettings=custom.settings,
                         gname = jgene))
  
  jstart <- subset(regions.annotated, region_coord == jpeak)$start - 6 * 10^5
  jend <- subset(regions.annotated, region_coord == jpeak)$end + 5 * 10^5
  
  
  # PlotGTrack(x.long %>% mutate(louvain = ifelse(louvain == 2, jtopic, "Others")),
  # jstart, jend, mart.obj, gen = "hg38", chr = "chr13", jheight = 1.5)
  
  # do hierarchical clustering on chromosome 13?
  # 
  # x.mat <- spread(x.long %>% dplyr::select(-c(louvain.orig, louvain)), key = cell, value = exprs)
  # rownames(x.mat) <- x.mat$coord
  # x.mat$coord <- NULL
  # clusters <- hclust(dist(x.mat))
  # plot(clusters)
  
  
  # Plot outputs to Dropbox -------------------------------------------------
  
  plotf <- paste0("~/Dropbox/scCHiC_figs/figures/K562_G1_Analysis/K562_G1_Analysis.mark.", jmark, ".", Sys.Date(), ".pdf")
  
  pdf(plotf, useDingbats = FALSE)
  
  jcounts <- Matrix::colSums(count.mat)
  jcol.counts <- ColorsByCounts(jcounts, nbreaks = 100, colvec = c("gray90", "gray50", "darkblue"))
  
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
  
  par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1))
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste(jmain, "color by total cell counts"), pty = 's', col = jcol.counts, asp = 0.5)
  
  par(mfrow=c(nb.row, nb.col), mar=c(1,0.5,0.5,1))
  mapply(function(jcol.rgb, jtopic){
    plot(dat.umap$layout[, 1], dat.umap$layout[, 2], pch = 20, main = paste("Topic", jtopic), col = jcol.rgb, asp = 0.5)
  }, jcol.rgbs, seq(kchoose))
  
  # plot umap with louvain
  m.louvain <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = as.character(louvain))) + geom_point() +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_brewer(palette = "Spectral")
  print(m.louvain)
  
  # plot without color
  m.nocolor <- ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(alpha = 0.5, size = 3) +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m.nocolor)
  
  PlotGTrack(x.long, jstart, jend, mart.obj, gen = "hg38", chr = jchromo, jheight = "auto")
  
  dev.off()
}

