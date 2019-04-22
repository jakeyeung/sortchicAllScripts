# Jake Yeung
# GetMetaData.R
#  
# 2018-12-19


GetTissue <- function(inf, type="BM"){
  # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> BM
  # PZ-K562-H3K27me3-G1_AH2VV5BGX9_S11.bc_counts.txt
  # PZ-K562-G1-M-H3K4me1-plate1_AHHCCGBGX9_S2 -> K562 ok
  if (type == "BM"){
    x <- strsplit(inf, split = "-")[[1]][[2]]
  } else if (type == "K562" | type == "K562_round2"){
    x <- strsplit(inf, split = "-")[[1]][[2]]
  } else {
    print("Must be BM or K562")
  return(x)
  }
}

GetChip <- function(inf, type="BM"){
  # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> H3K7me3
  # PZ-K562-H3K27me3-G1_AH2VV5BGX9_S11.bc_counts.txt -> H3K27me3
  # PZ-K562-G1-M-H3K4me1-plate1_AHHCCGBGX9_S2 -> H3K4me1
  if (type == "BM"){
    x <- strsplit(inf, split = "-")[[1]][[4]]
  } else if (type == "K562"){
    x <- strsplit(inf, split = "-")[[1]][[3]]
  } else if (type == "K562_round2"){
    x <- strsplit(inf, split = "-")[[1]][[5]]
  } else {
    print("Must be BM or K562")
  }
  return(x)
}

TechRepToRank <- function(cnames.new, mouse.reps = c("m1", "m2"), jsep = "-"){
  # expects m1 and m2 mouse replicates
  # sampname: PZ-BM-m1-H3K4me1-1_AH3VGVBGX9_S9
  for (mouse in mouse.reps){
    # get tech reps for each mouse
    cnames.i <- grep(mouse, cnames.new)
    if (length(cnames.i) == 0){
      warning(paste("Length of cnames after grepping mouse replicate is 0: Check the grep str", mouse))
    }
    cnames.sub <- cnames.new[cnames.i]
    techreps <- unique(sapply(cnames.sub, function(x) strsplit(x, split=jsep)[[1]][[4]]))  # repS9 and repS11
    techreps.number <- sapply(techreps, function(techrep) readr::parse_number(techrep))  # 9, 11
    # convert number to rank
    techreps.rank <- rank(techreps.number)
    # techreps.new <- paste("", techreps.rank, sep = "")  # can add prefix 
    techreps.new <- techreps.rank
    cnames.sub.new <- cnames.sub  # init
    for (i in seq(length(techreps.new))){
      techrep <- techreps[[i]]
      techrep.new <- techreps.new[[i]]
      cnames.sub.new <- gsub(paste0(jsep, techrep, jsep), paste0(jsep, techrep.new, jsep), cnames.sub.new)
    }
    # replace repS11 with rep1
    # cnames.sub.new <- mapply(function(techrep, techrep.new, cname) return(gsub(techrep, techrep.new, cname)), techreps, techreps.new, cnames.sub)
    cnames.new[cnames.i] <- cnames.sub.new
  }
  return(cnames.new)
}

GetBioRep <- function(inf, type="BM"){
  # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> m1
  # PZ-K562-H3K27me3-G1_AH2VV5BGX9_S11.bc_counts.txt -> G1
  # PZ-K562-G1-M-H3K4me1-plate1_AHHCCGBGX9_S2 -> plate1
  if (type == "BM"){
    # x <- strsplit(x, split = "m")[[1]][[2]]  # keep the m1
    x <- strsplit(inf, split = "-")[[1]][[3]]
  } else if (type == "K562"){
    xtmp <- strsplit(inf, split = "-")[[1]][[4]]
    x <- strsplit(xtmp, split = "_")[[1]][[1]]
  } else if (type == "K562_round2"){
    xtmp <- strsplit(inf, split = "-")[[1]][[6]]
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
  # if (type == "BM" | type == "K562" | type == "K562_round2"){
    x <- stringr::str_match(inf, "_S[0-9]*.")[1, 1]
    x <- substr(x, start = 2, stop = nchar(x) - 1)
  # }
  # if (type == "BM"){
  #   x <- stringr::str_match(inf, "-[1-9]*_")[1, 1] 
  #   # remove first and last characters
  #   x <- substr(x, start = 2, stop = nchar(x) - 1)
  # } else if (type == "K562"){
  #   x <- stringr::str_match(inf, "_S[1-9]*.")[1, 1]
  #   x <- substr(x, start = 2, stop = nchar(x) - 1)
  # }
  return(x)
}
GetExperiment <- function(inf, type="BM"){
  # "PZ-BM-m1-H3K27me3-1_H2GV2BGX9_S14.bc_counts.txt" -> H2GV2BGX9
  # PZ-K562-H3K27me3-G1_AH2VV5BGX9_S11.bc_counts.txt
  x <- strsplit(inf, split = "_")[[1]][[2]]
  return(x)
}
