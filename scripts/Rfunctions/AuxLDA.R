
DoLouvain <- function(topics.mat, custom.settings.louv, dat.umap.long = NULL){
  # Do Louvain for clustering
  dat.umap.louv <- umap(topics.mat, config = custom.settings.louv)
  dat.umap.louv.long <- data.frame(umap1 = dat.umap.louv$layout[, 1], umap2 = dat.umap.louv$layout[, 2], cell = rownames(dat.umap.louv$layout), 
                                   stringsAsFactors = FALSE)
  cell.indx <- hash(rownames(dat.umap.louv$knn$indexes), dat.umap.louv$knn$indexes[, 1])
  cell.indx.rev <- hash(dat.umap.louv$knn$indexes[, 1], rownames(dat.umap.louv$knn$indexes))
  nr <- nrow(dat.umap.louv$knn$indexes)
  nc <- ncol(dat.umap.louv$knn$indexes)
  edgelist <- matrix(NA, nrow = nr * nc, ncol = 2)
  colnames(edgelist) <- c("from", "to")
  for (vertex.i in seq(nr)){
    istart <- nc*(vertex.i - 1)+1
    iend <- nc*vertex.i
    edgelist[istart : iend, 1] <- cell.indx.rev[[as.character(vertex.i)]]
    edgelist[istart : iend, 2] <- sapply(dat.umap.louv$knn$indexes[vertex.i, 1:nc], function(x) cell.indx.rev[[as.character(x)]])
    # edgelist[istart : iend, 3] <- 1 / (dat.umap$knn$distances[vertex.i, 1:nc] + 0.1)
  }
  g <- graph_from_data_frame(edgelist, directed=FALSE)
  g.out <- cluster_louvain(g, weights = NULL)
  V(g)$color <- g.out$membership
  clstr <- hash(g.out$names, g.out$membership)
  if (is.data.frame(dat.umap.long)){
    dat.umap.long$louvain <- as.character(sapply(dat.umap.long$cell, function(x) clstr[[as.character(x)]]))
  } else {
    dat.umap.long <- clstr
  }
  return(dat.umap.long)
}

GetPeaksFromGene <- function(jgene, regions.annot, dist = 50000){
  jsub <- subset(regions.annot, grepl(jgene, SYMBOL)) %>% arrange(abs(distanceToTSS)) %>% filter(distanceToTSS <= dist)
  jpeaks <- jsub$region_coord
  return(list(regions.sub = jsub, peaks = jpeaks))
}

SelectBestPeak <- function(jpeaks, regions.annot = NULL, tm.result){
  if (length(jpeaks) == 1){
    warning("Only one peak entered, returning only option")
    # no need to select
    return(jpeaks)
  }
  terms.sub <- tm.result$terms[, jpeaks]
  jmax <- apply(terms.sub, 2, max)
  jmax.i <- apply(terms.sub, 2, which.max)  # which topic max occurs. Often they agree.
  
  print(unique(jmax.i))
  
  jpeak <- jpeaks[which(jmax == max(jmax))]
  return(jpeak)
}

PlotAllMarks <- function(jgene, jpeak, jmarks, out.objs, custom.settings){
  out <- lapply(jmarks, function(jmark){
    print(jmark)
    m <- PlotImputedPeaks(tm.result.lst[[jmark]], jpeak, jmark, show.plot = FALSE, return.plot.only = TRUE, usettings=custom.settings, gname = jgene)
  })
  return(out)
}


LoadLDABins <- function(jmark, jbin=TRUE, top.thres=0.995, inf = NULL, convert.chr20.21.to.X.Y = TRUE, add.chr.prefix = FALSE, choose.k = "auto"){
  # jbin <- "TRUE"
  # top.thres <- 0.995
  if (is.null(inf)){
    if (jbin){
      inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.", jbin, "/lda_out_meanfilt.BM-", jmark, ".CountThres0.K-5_10_15_20_25.Robj")
    } else {
      inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_MetaCell/lda_outputs.meanfilt_1.cellmin_100.cellmax_500000.binarize.", jbin, "/lda_out_meanfilt.BM-", jmark, ".CountThres0.K-5_15_25.Robj")
    }
  }
  assertthat::assert_that(file.exists(inf))
  load(inf, v=T)
  if (choose.k == "auto"){
    out.lda <- ChooseBestLDA(out.lda)
  } else {
    kvec <- lapply(out.lda, function(x) x@k)
    i <- which(kvec == as.numeric(choose.k))
    out.lda <- out.lda[[i]]
  }

  
  if (add.chr.prefix){
    out.lda@terms <- paste0("chr", out.lda@terms)
  }
  
  kchoose <- out.lda@k
  tm.result <- posterior(out.lda)
  
  if (convert.chr20.21.to.X.Y){
    colnames(tm.result$terms) <- gsub("chr20", "chrX", colnames(tm.result$terms))
    colnames(tm.result$terms) <- gsub("chr21", "chrY", colnames(tm.result$terms))
  }
  
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
  
  topic.regions <- lapply(seq(kchoose), function(clst){
    return(SelectTopRegions(tm.result$terms[clst, ], colnames(tm.result$terms), method = "thres", method.val = top.thres))
  })
  
  return(list('out.lda' = out.lda, 'tm.result' = tm.result, 'topic.regions' = topic.regions, 'regions.annotated' = regions.annotated, 'count.mat' = count.mat))
}

GetVar <- function(tm.result, regions.annotated){
  mat.norm <- t(tm.result$topics %*% tm.result$terms)
  # find peaks with largest range
  maxmin <- sort(apply(mat.norm, 1, function(jcol) mad(jcol)), decreasing = TRUE)
  # what are the sizes?
  peak.size <- sapply(names(maxmin), function(x) ParseCoord(x)$end - ParseCoord(x)$start)
  # normalize maxmin by peaksize
  maxmin.norm <- maxmin / peak.size
  dat.var <- data.frame(Var = maxmin, peak.size = peak.size, coord = names(maxmin))
  regions.annotated <- dplyr::left_join(regions.annotated, dat.var, by = c("region_coord"="coord"))
  return(regions.annotated)
}

ChooseBestLDA <- function(out.lda){
  # pick best k
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
  return(out.lda)
}


GetUmapSettings <- function(nn, jmetric, jmindist, seed=123){
  # nn <- 5
  # jmetric <- 'euclidean'
  # # jmetric <- 'cosine'
  # jmindist <- 0.1
  custom.settings <- umap.defaults
  custom.settings$n_neighbors <- nn
  custom.settings$metric <- jmetric
  custom.settings$min_dist <- jmindist
  custom.settings$random_state <- seed
  return(custom.settings)
}

# color by loadings on Kvec
ColorsByGamma <- function(topic, topics.mat, cols.vec = c("pink", "red", "darkred")){
  # jcol <- out.lda@gamma[, topic]
  # jcol <- tmResult$topics[, topic]
  jcol <- topics.mat[, topic]
  colorPal <- grDevices::colorRampPalette(cols.vec)
  jcol.rgb <- colorPal(200)[as.numeric(cut(jcol,breaks = 200))]
  return(jcol.rgb)
}

ColorsByCounts <- function(counts.vec, nbreaks=100, colvec = c("pink", "red", "darkred")){
  colorPal <- grDevices::colorRampPalette(colvec)
  jcol.rgb <- colorPal(nbreaks)[as.numeric(cut(counts.vec,breaks = nbreaks))]
  return(jcol.rgb)
}

SelectTopRegions <- function(beta.row, regions, method = "thres", method.val = 0.01){
  # take betas (in log scale) and select top X fraction
  # or do simple cutoff, set method = "cutoff"
  if (method == "cutoff"){
    return(regions[which(beta.row > method.val)])
  } else if (method == "thres"){
    return(regions[which(beta.row > quantile(beta.row, method.val))])
  } else {
    stop(paste("Method", method, "not yet implemented"))
  }
}

GetCountMatFromLDA <- function(out.lda){
  # https://stackoverflow.com/questions/20004493/convert-simple-triplet-matrixslam-to-sparse-matrixmatrix-in-r
  count.mat <- Matrix::sparseMatrix(i=out.lda@wordassignments$i, 
                                  j=out.lda@wordassignments$j, 
                                  x=out.lda@wordassignments$v, 
                                  dims=c(out.lda@wordassignments$nrow, out.lda@wordassignments$ncol))
  return(count.mat)
}


.modelMatSelection <- function(
  # from cisTopics package 
  object,
  target,
  method,
  all.regions=FALSE
){
  # Check info
  if (length(object@selected.model) < 1){
    stop('Please, run selectModel() first.')
  }
  
  if (target == 'cell'){
    if (method == 'Z-score'){
      modelMat <- scale(object@selected.model$document_expects, center=TRUE, scale=TRUE)
    }
    
    else if (method == 'Probability'){
      alpha <- object@calc.params[['runModels']]$alpha/length(object@selected.model$topic_sums)
      modelMat <- apply(object@selected.model$document_sums, 2, function(x) {(x + alpha)/sum(x + alpha)})
    }
    else{
      stop('Incorrect method selected. Chose method between "Z-score" and "Probability".')
    }
    colnames(modelMat) <- object@cell.names
    rownames(modelMat) <- paste0('Topic', 1:nrow(modelMat))
  }
  
  else if (target == 'region'){
    if (!all.regions){
      if (length(object@binarized.cisTopics) < 1){
        stop('Please, use binarizecisTopics() first for defining the high confidence regions for dimensionality reduction!')
      }
      else {
        regions <- unique(unlist(lapply(object@binarized.cisTopics, rownames)))
      }
    }
    
    topic.mat <- object@selected.model$topics
    
    if (method == 'NormTop'){
      normalizedTopics <- topic.mat/(rowSums(topic.mat) + 1e-05)
      modelMat <- apply(normalizedTopics, 2, function(x) x * (log(x + 1e-05) - sum(log(x + 1e-05))/length(x)))
    }
    
    else if (method == 'Z-score'){
      modelMat <- scale(object@selected.model$topics, center=TRUE, scale=TRUE)
    }
    
    else if (method == 'Probability'){
      beta <- object@calc.params[['runModels']]$beta
      topic.mat <- object@selected.model$topics
      modelMat <-  (topic.mat + beta)/rowSums(topic.mat + beta)
    }
    
    else{
      stop('Incorrect method selected. Chose "NormTop", "Z-score" and "Probability".')
    }
    
    colnames(modelMat) <- object@region.names
    rownames(modelMat) <- paste0('Topic', 1:nrow(modelMat))
    
    if (!all.regions){
      modelMat <- modelMat[,regions]
    }
  }
  
  else{
    stop('Please, provide target="cell" or "region".')
  }
  
  return(modelMat)
}
