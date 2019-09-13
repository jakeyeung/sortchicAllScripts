
fquad <- function(x, intercept, slope, curve){
  return(intercept + slope * x + curve * x^2)
}

GetMarkB6 <- function(x){
  return(strsplit(x, "-")[[1]][[4]])
}

GetRepB6 <- function(x){
  return(paste0("rep", strsplit(strsplit(x, "-")[[1]][[5]], "_")[[1]][[1]]))
}

GetCellB6 <- function(x, shift = 0){
  # shift by 1 to get cells to start at 1 rather than 0
  cell <- as.numeric(strsplit(x, "_")[[1]][[4]]) + shift
  return(paste0("cell", as.character(cell)))
}

GetCellB6.uniqname <- function(x, shift = 0){
  # shift by 1 to get cells to start at 1 rather than 0
  cell <- as.numeric(strsplit(x, "_")[[1]][[2]]) + shift
  return(paste0("cell", as.character(cell)))
}

ReadDinuc <- function(inf){
  dat <- fread(inf)[-1, ] %>%
    dplyr::rename(dinuc = sampleName) 
  dat.long <- gather(dat, key = "samp", value = "count", -dinuc) %>%
    rowwise() %>%
    mutate(count = ifelse(is.na(count), 0, count),
           mark = GetMarkB6(samp),
           repl = GetRepB6(samp),
           cell = GetCellB6.uniqname(samp, shift = 0))
  return(dat.long)
}

ReadDinucSplit <- function(inf){
  dat <- fread(inf)[-1, ] %>%
    dplyr::rename(dinuc = sampleName, strand = V2)
  dat.long <- gather(dat, key = "samp", value = "count", -dinuc, -strand) %>%
    rowwise() %>%
    mutate(count = ifelse(is.na(count), 0, count),
           mark = GetMarkB6(samp),
           repl = GetRepB6(samp),
           cell = GetCellB6(samp))
  return(dat.long)
}

ReadMat <- function(inf, as.sparse = TRUE){
  dat <- fread(inf)[-1, ] %>%
    dplyr::rename(chromo = sampleName,
                  startend = V2)
  # remove missing
  dat <- subset(dat, startend != "Missing")
  dat$coord <- paste(dat$chromo, dat$startend, sep = ":")
  dat$chromo <- NULL
  dat$startend <- NULL
  if (as.sparse){
    coords <- dat$coord
    dat$coord <- NULL
    dat <- as.matrix(dat)
    dat[is.na(dat)] <- 0
    rownames(dat) <- coords
    dat <- Matrix::Matrix(dat, sparse = TRUE)
  }
  return(dat)
}

ReadCellSum <- function(inf){
  dat <- fread(inf)
}



# Linneg ------------------------------------------------------------------

GetMarkLinneg <- function(x){
  return(strsplit(x, "-")[[1]][[5]])
}

GetRepLinneg <- function(x){
  return(paste0("rep", strsplit(strsplit(x, "-")[[1]][[6]], "_")[[1]][[1]]))
}

GetCellLinneg <- function(x, shift = 0){
  cell <- as.numeric(strsplit(x, "_")[[1]][[4]]) + shift
  assertthat::assert_that(!is.na(cell))
  return(paste0("cell", as.character(cell)))
}

ReadDinucLinneg <- function(inf){
  dat <- fread(inf)[-1, ] %>%
    dplyr::rename(dinuc = sampleName) 
  dat.long <- gather(dat, key = "samp", value = "count", -dinuc) %>%
    rowwise() %>%
    mutate(count = ifelse(is.na(count), 0, count),
           mark = GetMarkLinneg(samp),
           repl = GetRepLinneg(samp),
           cell = GetCellLinneg(samp, shift = 0))
  return(dat.long)
}



# Get empty wells ---------------------------------------------------------


# #  1 index
GetEmptyWells <- function(indx = 0){
  if (indx == 1){
    indx.all <- seq(384)
    hascell.indx <- c(seq(1:356),seq(360:379)+360)
    empty.indx <- setdiff(indx.all, hascell.indx)
    empty.names <- paste("cell", empty.indx, sep = "")
  } else if (indx == 0){
    # 0 index
    indx.all <- seq(384) - 1
    hascell.indx <- c(seq(1:356),seq(360:379)+360) - 1
    empty.indx <- setdiff(indx.all, hascell.indx)
    empty.names <- paste("cell", empty.indx, sep = "")
  } else {
    warning("Index must be 0 or 1")
  }
  return(empty.names)
}


# For loading matrices ----------------------------------------------------


LoadMats <- function(infs, cells.keep){
  sparse.mats <- lapply(infs, function(inf) readRDS(inf))
  # filter cells
  # sparse.mats <- lapply(sparse.mats, function(mat){
  #   cnames <- colnames(mat)
  #   cnames.keep.i <- which(cnames %in% cells.keep)
  #   return(mat[, cnames.keep.i])
  # })
  rnames <- lapply(sparse.mats, rownames)
  rnames.common <- Reduce(intersect, rnames)
  sparse.mats <- lapply(sparse.mats, function(x){
    return(x[rnames.common, ])
  })
  mat.merge <- do.call(cbind, sparse.mats)
  cells.keep.i <- which(colnames(mat.merge) %in% cells.keep)
  mat.merge <- mat.merge[, cells.keep.i]
  return(mat.merge)
}

LoadBlacklist <- function(inf = "data/blacklists/mm10.blacklist.bed.gz", asGR = TRUE){
  dat <- fread(inf, col.names = c("seqnames", "start", "end"))
  if (asGR){
    dat <- makeGRangesFromDataFrame(dat)
  }
  return(dat)
}


