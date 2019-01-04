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
  xbin <- Matrix::Matrix(x, sparse = TRUE, nrow = nrow(x), ncol = ncol(x))
  return(xbin)
}
