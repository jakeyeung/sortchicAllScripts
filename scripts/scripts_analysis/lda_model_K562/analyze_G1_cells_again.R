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
colsvec <- list(H3K4me1 = "cyan1", H3K4me3 = "darkblue", H3K9me3 = "red1", H3K27me3 = "darkorange1")

for (jmark in jmarks){
  
  jcol <- colsvec[[jmark]]
  print(jcol)
  
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
  # 
  # top.thres <- 0.995
  # topic.regions <- lapply(seq(kchoose), function(clst){
  #   return(SelectTopRegions(tm.result$terms[clst, ], colnames(tm.result$terms), method = "thres", method.val = top.thres))
  # })
  # top.regions <- unique(unlist(topic.regions))
  # 
  # # Process regions ---------------------------------------------------------
  # 
  # # chainpath <- "~/data/scchic/databases/hg38ToHg19.over.chain"
  # # assertthat::assert_that(file.exists(chainpath))
  # # ch <- rtracklayer::import.chain(chainpath)
  # # seqlevelsStyle(regions.range) = "UCSC"
  # 
  # regions <- data.frame(seqnames = sapply(colnames(tm.result$terms), GetChromo),
  #                       start = sapply(colnames(tm.result$terms), GetStart),
  #                       end = sapply(colnames(tm.result$terms), GetEnd),
  #                       terms = colnames(tm.result$terms), 
  #                       stringsAsFactors = FALSE) %>%
  #   filter(seqnames != "chr23" &  start != "")  # negative values are returned as "", remove them
  # 
  # 
  # rownames(regions) <- regions$terms
  # 
  # regions.range <- makeGRangesFromDataFrame(as.data.frame(regions)) 
  # 
  # regions.annotated <- as.data.frame(annotatePeak(regions.range,
  #                                                 TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene,
  #                                                 annoDb='org.Hs.eg.db'))
  # regions.annotated <- regions.annotated %>%
  #   rowwise() %>%
  #   mutate(region_coord = paste(paste(seqnames, start, sep = ":"), end, sep = "-"))
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
  
  jtopic <- 1
  topn <- 50
  # find top hits
  jsub <- subset(top.peaks, topic == jtopic)
  
  jpeaks <- jsub$term[1:topn]
  
  print(jpeaks)
    
  print(head(subset(regions.annotated, region_coord %in% jpeaks)))
  
  # Do louvain on this ------------------------------------------------------
  
  dat.umap.long <- data.frame(umap1 = dat.umap$layout[, 1], umap2 = dat.umap$layout[, 2], cell = rownames(dat.umap$layout), stringsAsFactors = FALSE)
  
  # plot without color
  
  # Plot outputs to Dropbox -------------------------------------------------
  
  plotf <- paste0("~/Dropbox/scCHiC_figs/figures/K562_G1_Analysis/K562_G1_Analysis.mark.", jmark, ".", Sys.Date(), ".pdf")
  
  pdf(plotf, useDingbats = FALSE)
    m.nocolor <- ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point(alpha = 0.5, size = 3, color = jcol) +
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            axis.ticks=element_blank(), 
            axis.text.x=element_blank(),
            axis.text.y=element_blank(),
            panel.border=element_blank()) + 
      xlab("") + ylab("")
    print(m.nocolor)
  dev.off()
}


