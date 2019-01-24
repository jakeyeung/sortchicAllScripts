# Jake Yeung
# Date of Creation: 2019-01-23
# File: ~/projects/scChiC/scripts/scripts_analysis/lda_model/summarize_histone_marks_compare_metacell_binned.R
# Binned analysis 

rm(list=ls())

library(topicmodels)
library(dplyr)
library(ggplot2)
library(umap)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(rGREAT)
library(hash)

source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")

GetBins <- function(jchromo, midpt, jname, Q, winsize = 100000L){
  jleft <- midpt - winsize / 2
  jright <- midpt + winsize / 2
  Qsub <- subset(Q, chromo == jchromo & coord > jleft & coord < jright)
  Qsub$name <- jname
  return(as.data.frame(Qsub))
}

Vectorize(AssignHash <- function(x, jhash, null.fill = NA){
  # assign hash key to hash value, handle NULLs
  # null.fill = "original", returns original value x into jhash
  x.mapped <- jhash[[as.character(x)]]
  if (is.null(x.mapped)){
    if (is.na(null.fill)){
      x.mapped <- null.fill
    } else if (as.character(null.fill) == "original"){
      x.mapped <- x
    } else {
      x.mapped <- null.fill
    }
  }
  return(x.mapped)
}, vectorize.args = "x")

GetGOData <- function(out.tb.lst, ontology, jdecreasing=TRUE, order.by="Hyper_Fold_Enrichment"){
  topics <- which(lapply(1:length(out.tb.lst), function(i) !is.null(out.tb.lst[[i]][[ontology]])) == TRUE)
  GOdata <- lapply(topics, function(i) out.tb.lst[[i]][[ontology]])
  GOdata <- lapply(1:length(GOdata), function(i) GOdata[[i]][order(GOdata[[i]][order.by], 
                                                                   decreasing = jdecreasing),])
  GOdata <- lapply(1:length(GOdata), function(i) GOdata[[i]][1:top,])
  
  # annotate each GOdata by topic number
  for (i in topics){
    GOdata[[i]]$topic <- i
  }
  GOdata <- dplyr::bind_rows(GOdata) %>%
    mutate(topic = factor(topic, levels = sort(unique(topic))),
           onto = ontology)
  return(GOdata)
}



# Parameters to get files -------------------------------------------------


outdir <- "~/Dropbox/scCHiC_figs/FIG4_BM/LDA_outputs/binned"
dir.create(outdir)


jchip <- "H3K9me3"
jchip <- "H3K27me3"
jchip <- "H3K4me1"
jchip <- "H3K4me3"
  
  # jdist <- 1000L
  # jmean <- 1
  # jmin <- 100L
  # jmax <- 500000L
  # binarize <- "TRUE";  jtops <- "5_7_10_12_15_20_25_30"
  # binarize <- "FALSE"; jtops <- "15_20_25_30_35"
  # binarize <- "TRUE"; jtops <- "15_20_25_30_35"
  
  # jdir <- paste0('/tmp/ldaAnalysisHiddenDomains_', jdist, '/lda_outputs.meanfilt_', jmean, '.cellmin_', jmin, '.cellmax_', jmax, '.binarize.', binarize)
  # inf <- file.path(jdir, paste0('lda_out_meanfilt.PZ-BM-', jchip, '.CountThres0.K-', jtops, '.Robj'))
  inf <- paste0("/private/tmp/lda_output_binned/lda_out_meanfilt.BM-", jchip, ".CountThres0.K-5_10_15_20_25.Robj")
  inf.mc <- file.path(paste0("~/Dropbox/scCHiC_figs/FIG4_BM/", jchip, ".datadir_mc_f.Rda"))
  
  assertthat::assert_that(file.exists(inf))
  assertthat::assert_that(file.exists(inf.mc))
  
  load(inf, v=T)
  
  
  # what's the sparsity?
  print(1 - Matrix::nnzero(count.mat) / length(count.mat))
  

  
  # if already loaded the GREAT object, then don't need to choose best K
  
  if (length(out.lda) > 1){
    out.lda.lst <- out.lda
    # we did multicore, so which one to choose?
    Kvec <- sapply(out.lda.lst, function(x) x@k)
    best.K <- Kvec[which.max(sapply(out.lda.lst, function(x) x@loglikelihood))]
    # plot loglikelihood
    par(mfrow=c(1, 1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
    plot(Kvec, sapply(out.lda.lst, function(x) x@loglikelihood), 'o')
    kchoose <- best.K
    out.lda <- out.lda.lst[[which(Kvec == kchoose)]]
    print(paste("Likelihood: ", out.lda@loglikelihood))
  } else {
    kchoose <- out.lda@k
  }
  
  tm.result <- posterior(out.lda)
  top.thres <- 0.96 
  topic.regions <- lapply(seq(kchoose), function(clst){
    return(SelectTopRegions(tm.result$terms[clst, ], colnames(tm.result$terms), method = "thres", method.val = top.thres))
  })
  # regions <- data.frame(seqnames = sapply(colnames(tm.result$terms), GetChromo),
  #                       start = sapply(colnames(tm.result$terms), GetStart),
  #                       end = sapply(colnames(tm.result$terms), GetEnd))
  # regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
  # regions.annotated <- as.data.frame(annotatePeak(regions.range, 
  #                                                 TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene, 
  #                                                 annoDb='org.Mm.eg.db'))
  # regions.annotated$region_coord <- rownames(regions.range)

  
  
  # print(head(tm.result$topics))
  # cluster the terms, but first get interesting terms
  top.regions <- unique(unlist(topic.regions))  # topic.regions defined by threshold, 98 or 96th percentile of top weights in each column of the betas matrix
  
  terms.mat <- t(tm.result$terms)[top.regions, ]
  
  # Compare with metacell ---------------------------------------------------
  
  
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
  
  
  # Plotting everything now ----------------------------------------------------------------
  

  # # UMAP settings 
  # nn <- 5
  # jmetric <- 'euclidean'
  # # jmetric <- 'cosine'
  # jmindist <- 0.1
  # custom.settings <- umap.defaults
  # custom.settings$n_neighbors <- nn
  # custom.settings$metric <- jmetric
  # custom.settings$min_dist <- jmindist

  nn=5
  nnterms=15
  jmetric='euclidean' 
  jmindist=0.25
  custom.settings <- GetUmapSettings(nn=nn, jmetric=jmetric, jmindist=jmindist)
  custom.settings.terms <- GetUmapSettings(nn=nnterms, jmetric=jmetric, jmindist=jmindist)
  # # use different settings for peak weights for H3K4me1 and H3K4me3 for some reason
  # if (jchip %in% c("H3K4me1", "H3K4me3")){
  #   custom.settings.terms <- GetUmapSettings(nn=15, jmetric='euclidean', jmindist=0.01)
  # } else {
  #   custom.settings.terms <- custom.settings
  # }

  dat.umap <- umap(topics.mat, config = custom.settings)
  rownames(dat.umap$layout) <- rownames(topics.mat)
  jmain <- paste("Neighbors", nn, "Metric", jmetric, "MinDist", jmindist)
  
  # # run umap again on all cells or just common between MetaCell and LDA?
  # dat.umap.all <- umap(tm.result$topics, config = custom.settings)
  # rownames(dat.umap.all$layout) <- rownames(tm.result$topics)
  
  # compare metacell and lda
  plot(dat.umap$layout[, 1], dat.umap$layout[, 2], 
       pch = 20, 
       main = paste0("MC colors: ", jchip), 
       col = cells.rgb, 
       asp = 1, 
       xlab = "UMAP Dim 1", ylab = "UMAP Dim 2")
  