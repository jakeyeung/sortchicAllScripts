# Jake Yeung
# Date of Creation: 2019-05-13
# File: ~/projects/scchic/scripts/Rfunctions/AuxB6.R
# Some B6-specific functions 


LoadGetImputed <- function(inf){
  load(inf)  # out.objs, dat.umap.long, custom.settings. louv.settings
  imput.mat <- t(out.objs$tm.result$topics %*% out.objs$tm.result$terms)
  return(imput.mat)
}

LoadGetTmResult <- function(inf){
  load(inf)  # out.objs, dat.umap.long, custom.settings. louv.settings
  return(out.objs$tm.result)
}

LoadUmap <- function(inf.dat){
  assertthat::assert_that(file.exists(inf.dat))
  print(paste("Loading from:", inf.dat))
  load(inf.dat, v=F)
  dat.umap.long$repl <- as.character(sapply(as.character(dat.umap.long$cell), function(x) strsplit(strsplit(x, "-")[[1]][[5]], "_")[[1]][[1]]))
  dat.umap.long$techname <- as.character(sapply(as.character(dat.umap.long$cell), function(x) strsplit(strsplit(x, "-")[[1]][[5]], "_")[[1]][[1]]))
  dat.umap.long$cellindx <- as.character(sapply(as.character(dat.umap.long$cell), function(x) strsplit(strsplit(x, "-")[[1]][[5]], "_")[[1]][[2]]))
  return(dat.umap.long)
}


LoadFACS <- function(jmark, repl, cnames.noX, prefix){
  inf.facs <- paste0("/Users/yeung/data/scchic/facs/20190403/20190403_Mouse_BM_B6_", jmark, "_00", repl, "_index.csv")
  dat.facs <- read.csv(inf.facs)
  dat.facs$X <- NULL  # confuses downstream colnames
  dat.facs$repl <- as.character(repl)
  # dat.facs$techname <- sapply(as.character(dat.facs$repl), function(x) rep2techname[[x]])
  dat.facs$mark <- jmark
  dat.facs$wellindx <- as.numeric(c(seq(1:356), seq(360:379)+360))
  dat.facs$wellindx0 <- dat.facs$wellindx - 1  # zero based
  dat.facs$cell <- paste0(prefix, dat.facs$repl, "_", dat.facs$wellindx0)
  # keep colnames that start with X or in cnamesnoX
  cnames.keep <- colnames(dat.facs)[grepl(pattern = "^X", colnames(dat.facs))]
  # add noX
  cnames.keep <- c(cnames.noX, cnames.keep)
  dat.facs.filt <- dat.facs %>% dplyr::select(cell, cnames.keep)
  return(dat.facs.filt)
}

GetRefCnames <- function(inf.facs.ref = "/Users/yeung/data/scchic/facs/H3K9me3_index_2mice_4plates.csv"){
  # inf.facs.ref <- 
  dat.refs <- read.table(inf.facs.ref)
  cnames.noX <- colnames(dat.refs)[!grepl(pattern = "^X", colnames(dat.refs))]
}

GetFACSloadings <- function(dat.facs.filt, PC=1, make.pos = FALSE){
  X <- prcomp(dat.facs.filt, center = TRUE, scale. = TRUE)
  xout <- X$x[, PC]
  if (make.pos){
    if (xout[which.max(abs(xout))] < 0){
      xout <- xout * -1
    }
  }
  return(xout)
}

LoadFACSGetLoadings <- function(jmark, prefix){
  if (missing(prefix)){
    prefix <- paste0("B6-13W1-BM-", jmark, "-")
  }
  if (jmark != "H3K4me1"){
    repl.vec <- seq(4)
  } else {
    repl.vec <- c(2, 3, 4)
  }
  cnames.noX <- GetRefCnames()
  dat.facs.filt <- lapply(repl.vec, function(repl) LoadFACS(jmark, repl, cnames.noX, prefix)) %>%
    bind_rows()
  rownames(dat.facs.filt) <- dat.facs.filt$cellname
  dat.pca <- GetFACSloadings(dat.facs.filt %>% dplyr::select(-cell), PC = seq(5))
  dat.facs.filt$loadings <- GetFACSloadings(dat.facs.filt %>% dplyr::select(-cell), PC = 1, make.pos = TRUE)
  dat.facs.filt$mark = jmark
  return(dat.facs.filt)
}