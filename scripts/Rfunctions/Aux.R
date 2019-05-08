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


GetPeakSize <- function(coord){
  # chr1:3005258-3006803 -> 1545
  jstart <- as.numeric(strsplit(strsplit(coord, ":")[[1]][[2]], "-")[[1]][[1]])
  jend <- as.numeric(strsplit(strsplit(coord, ":")[[1]][[2]], "-")[[1]][[2]])
  return(jend - jstart)
}

CleanCoords <- function(x){
  # "chr11:3,089,532-3,115,355" -> remove commas
  return(gsub(pattern = ",", "", x))
}

BinarizeMatrix <- function(x){
  # https://stackoverflow.com/questions/14526429/turn-a-count-matrix-into-a-binary-existence-matrix
  xbin <- as.numeric(as.matrix(x) > 0)
  xbin <- Matrix::Matrix(xbin, sparse = TRUE, nrow = nrow(x), ncol = ncol(x))
  # get back the column and row names
  rownames(xbin) <- rownames(x)
  colnames(xbin) <- colnames(x)
  return(xbin)
}

ParseCoord <- function(x){
  # chr7:103,796,583-103,857,605 -> chromo, start, end
  out <- list()
  out$chromo <- strsplit(x, split = ":")[[1]][[1]]
  out$start <- strsplit(strsplit(x, split = ":")[[1]][[2]], split = "-")[[1]][[1]]
  out$end <- strsplit(strsplit(x, split = ":")[[1]][[2]], split = "-")[[1]][[2]]
  # remove commas
  out <- lapply(out, function(x) gsub(",", "", x))
  out$start <- as.numeric(gsub(",", "", out$start))
  out$end <- as.numeric(gsub(",", "", out$end))
  return(out)
}

GetChromo <- function(x, add.chr = FALSE){
  # chrY:90799295-90803056 -> chrY
  chromo <- strsplit(x, ":")[[1]][[1]]
  if (add.chr){
    chromo <- paste0("chr", chromo)
  }
  return(chromo)
}

GetStart <- function(x){
  # chrY:90799295-90803056 -> 90799295
  return(strsplit(strsplit(x, ":")[[1]][[2]], "-")[[1]][[1]])
}

GetEnd <- function(x){
  # chrY:90799295-90803056 -> 90799295
  return(strsplit(strsplit(x, ":")[[1]][[2]], "-")[[1]][[2]])
}
