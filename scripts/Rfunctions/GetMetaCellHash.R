# Jake Yeung
# Date of Creation: 2019-01-09
# File: ~/projects/scChiC/scripts/Rfunctions/GetMetaCellHash.R
# Functions

GetReplHash <- function(experinames.m, biorep, jhash){
  # get replicate names from each mouse
  gstr <- paste0('-BM-', biorep)  # add m?
  samps <- grep(gstr, experinames.m, value = TRUE)
  replval.old <- as.numeric(sapply(samps, function(x){
    xtmp <- strsplit(x, "-")[[1]][[4]]
    return(gsub("S", "", xtmp))
  }))
  replval.new <- paste("rep", rank(replval.old), sep = "")
  for (i in seq(length(samps))){
    jhash[samps[i]] <- replval.new[i]
  }
  return(jhash)
}

MakeNewExperiName <- function(old.name, replhash){
  repval <- replhash[[old.name]]
  xshort <- paste(strsplit(old.name, "-")[[1]][c(2,1,3)], collapse = "-")  # switch BM with H3K27me3
  xnew <- paste0(xshort, "-", repval)
  return(xnew)
}

  MakeNewCellName.rev <- function(x, experihash, cellhash, return.path=TRUE){
    # BM_H3K4me1_m1_rep1_cell4 -> H3K27me3-BM-1-S14-H2GV2BGX9-AACTAGTG
    # goal: match file names here: /hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/bam_split_by_bc/count_thres-0/PZ-BM-m1-H3K4me1-1_AH3VGVBGX9_S9
    indx <- gsub("cell", "", strsplit(x, "_")[[1]][[5]])
    bc <- AssignHash(indx, cellhash)
    jsci <- "PZ"
    jtiss <- strsplit(x, "_")[[1]][[1]]
    jmark <- strsplit(x, "_")[[1]][[2]]
    jmouse <- strsplit(x, "_")[[1]][[3]]
    jrepl <- gsub("rep", "", strsplit(x, "_")[[1]][[4]])
    jkey <- paste(jsci, jtiss, jmouse, jmark, jrepl, sep = "-")
    jval <- AssignHash(jkey, experihash)
    assertthat::assert_that(!is.na(jval))
    # link key and val to make dirname
    dirname <- paste(jkey, jval, sep = "_")
    # make filename
    # "PZ-BM-m1-H3K4me1-1_AH3VGVBGX9_S9.filtered.sorted.CTCTGAAG.sorted.bam.bai"
    fname <- paste0(dirname, ".filtered.sorted.", bc, ".sorted.bam")
	if (return.path){
	  out <- file.path(dirname, fname)
	} else {
	  out <- fname
	}
	return(out)
  }

MakeNewCellName <- function(oldname, replhash, cellhash, add.m = TRUE){
  # H3K27me3-BM-1-S14-H2GV2BGX9-AACTAGTG -> BM-H3K27me3-m1-rep2-cell244
  # replhash: convert old experiment name to new experiment name
  #
  
  # add m to biological replicate
  if (add.m){
    oldname <- gsub("-BM-", "-BM-m", oldname)
  }
  experi <- paste(strsplit(oldname, "-")[[1]][-6], collapse = "-")
  cellname <- strsplit(oldname, "-")[[1]][[6]]
  experi.new <- MakeNewExperiName(experi, replhash)
  cellname.new <- cellhash[[cellname]]
  newname <- paste0(experi.new, "-", cellname.new)
  # replace dashes with underscores
  newname <- gsub("-", "_", newname)
  return(newname)
}

GetMCCoords <- function(object, winsize = 100000){
  Q<-as.data.frame(object@mc_fp)
  rownames(Q)<-paste('chr',rownames(Q),sep='')  # easier to grep with chromosome as start
  Q$fullcoord <- rownames(Q)
  Q <- data.frame(fullcoord = as.character(Q$fullcoord), stringsAsFactors = FALSE)
  Q <- Q %>%
    rowwise() %>%
    mutate(coord = as.integer(strsplit(fullcoord, split = "_")[[1]][[2]]) + (winsize/2),
           chromo = strsplit(fullcoord, split = "_")[[1]][[1]],
           left = coord - winsize / 2,
           right = coord + winsize / 2)
  return(Q)
}

GetMCObjects <- function(inf){
  load(inf, v=F)  # object
  # Use same coloring as MetaCell output so we can directly compare
  Q <- GetMCCoords(object)
  mc_index<-object@mc
  # mc_colors <- sapply(object@mc, function(x) ifelse(x %in% c(5, 7), "#A9A9A9", "#000000"))
  mc_colors<-object@colors
  return(list("mc_index" = mc_index, "mc_colors" = mc_colors, "Q" = Q, "exprs" = object@mc_fp))
}

GetQLong <- function(inf){
  load(inf, v=F)
  Q<-as.data.frame(object@mc_fp)
  rownames(Q)<-paste('chr',rownames(Q),sep='')  # easier to grep with chromosome as start
  Q$fullcoord <- rownames(Q)
  Qlong <- reshape2::melt(Q, value.name = "exprs", variable.name = "MC", id.vars = "fullcoord")
  colnames(Qlong) <- c("fullcoord", "MC", "exprs")
  Qlong$MC <- paste0("mc", Qlong$MC)
  return(Qlong)
}
