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