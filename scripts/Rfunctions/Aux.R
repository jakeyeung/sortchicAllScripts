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

GetChromo <- function(x){
  # chrY:90799295-90803056 -> chrY
  return(strsplit(x, ":")[[1]][[1]])
}

GetStart <- function(x){
  # chrY:90799295-90803056 -> 90799295
  return(strsplit(strsplit(x, ":")[[1]][[2]], "-")[[1]][[1]])
}

GetEnd <- function(x){
  # chrY:90799295-90803056 -> 90799295
  return(strsplit(strsplit(x, ":")[[1]][[2]], "-")[[1]][[2]])
}
