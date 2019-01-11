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
