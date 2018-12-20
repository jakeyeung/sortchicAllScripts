# Jake Yeung
# GetMetaData.R
#  
# 2018-12-19

GetTissue <- function(inf){
  # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> BM
  x <- strsplit(inf, split = "-")[[1]][[2]]
  return(x)
}
GetChip <- function(inf){
  # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> H3K7me3
  x <- strsplit(inf, split = "-")[[1]][[4]]
  return(x)
}
GetBioRep <- function(inf){
  # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> m1 -> 1
  x <- strsplit(inf, split = "-")[[1]][[3]]
  x <- strsplit(x, split = "m")[[1]][[2]]
  return(x)
}
GetTechRep <- function(inf){
  # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> S14-> 14
  x <- stringr::str_match(inf, "_S[1-9]*.")[1, 1]
  # remove first and last characters
  x <- substr(x, start = 2, stop = nchar(x) - 1)
  return(x)
}
GetExperiment <- function(inf){
  # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> H2GV2BGX9
  x <- strsplit(inf, split = "_")[[1]][[2]]
  return(x)
}
