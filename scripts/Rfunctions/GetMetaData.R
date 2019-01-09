# Jake Yeung
# GetMetaData.R
#  
# 2018-12-19

GetTissue <- function(inf, type="BM"){
  # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> BM
  # PZ-K562-H3K27me3-G1_AH2VV5BGX9_S11.bc_counts.txt
  if (type == "BM"){
    x <- strsplit(inf, split = "-")[[1]][[2]]
  } else if (type == "K562"){
    x <- strsplit(inf, split = "-")[[1]][[2]]
  } else {
    print("Must be BM or K562")
  }
  return(x)
}
GetChip <- function(inf, type="BM"){
  # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> H3K7me3
  # PZ-K562-H3K27me3-G1_AH2VV5BGX9_S11.bc_counts.txt
  if (type == "BM"){
    x <- strsplit(inf, split = "-")[[1]][[4]]
  } else if (type == "K562"){
    x <- strsplit(inf, split = "-")[[1]][[3]]
  } else {
    print("Must be BM or K562")
  }
  return(x)
}
GetBioRep <- function(inf, type="BM"){
  # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> m1
  # PZ-K562-H3K27me3-G1_AH2VV5BGX9_S11.bc_counts.txt -> G1
  if (type == "BM"){
    # x <- strsplit(x, split = "m")[[1]][[2]]  # keep the m1
    x <- strsplit(inf, split = "-")[[1]][[3]]
  } else if (type == "K562"){
    xtmp <- strsplit(inf, split = "-")[[1]][[4]]
    x <- strsplit(xtmp, split = "_")[[1]][[1]]
  } else {
    print("Must be BM or K562")
  }
  return(x)
}
GetTechRep <- function(inf, type="BM"){
  # # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> S14-> 14
  # PZ-K562-H3K27me3-G1_AH2VV5BGX9_S11.bc_counts.txt
  # x <- stringr::str_match(inf, "_S[1-9]*.")[1, 1]
  # "PZ-BM-m2-H3K9me3-1_H2GV2BGX9_S17" -> 1
  if (type == "BM"){
    x <- stringr::str_match(inf, "-[1-9]*_")[1, 1] 
    # remove first and last characters
    x <- substr(x, start = 2, stop = nchar(x) - 1)
  } else if (type == "K562"){
    x <- stringr::str_match(inf, "_S[1-9]*.")[1, 1]
    x <- substr(x, start = 2, stop = nchar(x) - 1)
  }
  return(x)
}
GetExperiment <- function(inf, type="BM"){
  # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> H2GV2BGX9
  # PZ-K562-H3K27me3-G1_AH2VV5BGX9_S11.bc_counts.txt
  x <- strsplit(inf, split = "_")[[1]][[2]]
  return(x)
}
