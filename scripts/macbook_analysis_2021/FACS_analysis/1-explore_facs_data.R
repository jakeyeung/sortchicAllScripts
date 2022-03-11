# Jake Yeung
# Date of Creation: 2022-01-25
# File: ~/projects/scchic/scripts/macbook_analysis_2021/FACS_analysis/1-explore_facs_data.R
# 
# copy from 
# File: ~/projects/dblchic/scripts/macbook_analysiis/integrate_with_FACS/get_umap_in_FACS_space.R

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(topicmodels)

library(hash)
library(igraph)
library(umap)

library(ggrastr)


# Functions ---------------------------------------------------------------


GetAnnotFromName <- function(fullname){
  # remove all "04-01 09-03"
  # batch <- strsplit(basename(fullname), split = " ")[[1]][[3]]
  is.dbl <- stringi::stri_detect(str = fullname, regex = " 04-01 09-03 ")
  
  batch <- stringi::stri_extract(str = fullname, regex = "SL[1-9]")
  if (is.dbl){
    mark <- "k4me1-k9me3"
  } else {
    mark <- stringi::stri_extract(fullname, regex = "\\d+-\\d+")
  }
  platelong <- stringi::stri_extract(fullname, regex = "\\d+_index.csv$")
  plate <- strsplit(platelong, split = "_")[[1]][[1]]
  # if (batch == "SL3"){
  #   # BM SL3 04-01 Plate_001_index.csv
  #   mark <- strsplit(basename(fullname), split = " ")[[1]][[4]]
  #   platelong <- strsplit(basename(fullname), split = " ")[[1]][[5]]
  #   assertthat::assert_that(startsWith(platelong, "Plate"))
  #   plate <- strsplit(platelong, split = "_")[[1]][[2]]
  # } else if (batch == "SL4"){
  #   # 20211125 Mouse BM SL4 4-3_001_index.csv
  #   mark <- strsplit(basename(fullname), split = " ")[[1]][[4]]
  #   platelong <- strsplit(basename(fullname), split = " ")[[1]][[5]]
  #   assertthat::assert_that(startsWith(platelong, "Plate"))
  #   plate <- strsplit(platelong, split = "_")[[1]][[2]]
  # }
  return(list(batch = batch, mark = mark, plate = plate))
}



GetMarkToCode <- function(){
  jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
  jcodes <- c("04-01", "04-03", "27-03", "09-03"); names(jcodes) <- jcodes
  
  jmarks2 <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
  jcodes2 <- c("04-01", "04-03", "27-03", "09-03"); names(jcodes) <- jcodes
  
  hash::hash(jmarks, jcodes)
}

WellToCellNumber <- function(w, totalcols = 24, shift = -1){
  # Well name (e.g., P20 to Plate coordinate (e.g. 380)
  jlet <- toupper(letters[1:26])
  jrow <- as.numeric(match(gsub("[^a-zA-Z]", "", w), jlet))
  jcol <- as.numeric(gsub("[a-zA-Z]", "", w))
  # shift -1 gives zerobased coordinates
  return( as.character((jrow - 1) * totalcols + jcol + shift) )
}

WellToRow <- function(w){
  # Well name (e.g., P20 to Plate coordinate (e.g. 380)
  jlet <- toupper(letters[1:26])
  jrow <- as.numeric(match(gsub("[^a-zA-Z]", "", w), jlet))
  return(jrow)
}

WellToColumn <- function(w){
  # Well name (e.g., P20 to Plate coordinate (e.g. 380)
  jlet <- toupper(letters[1:26])
  jcol <- as.numeric(gsub("[a-zA-Z]", "", w))
  return(jcol)
}

GetFiles <- function(indir, mark = "k4me1", batch = "SL3"){
  # BM SL3 04-01 Plate_001_index.csv
  m2c <- GetMarkToCode()
  code <- m2c[[mark]]
  jpat <- paste0(" ", batch, " ", code, " Plate_.*.csv")
  print(jpat)
  list.files(path = indir, all.files = TRUE, full.names = TRUE, pattern = jpat)
}

GetFiles2 <- function(indir, batch = "SL3"){
    # BM SL3 04-01 Plate_001_index.csv
    # jpat <- paste0(" ", batch, ".*.Plate_.*.csv")
    jpat <- paste0(" ", batch, ".*.csv")
    print(jpat)
    infs <- list.files(path = indir, all.files = TRUE, full.names = TRUE, pattern = jpat)
  return(infs)
}

DoUmapOnFACS <- function(indir, jplate.i){
  jpat <- paste0(".*All K.*plate ", jplate.i, ".*.csv")
  print(jpat)
  inf.facspath <- list.files(indir, pattern = jpat, full.names = TRUE)
  print(inf.facspath)
  inf.facsname <- basename(inf.facspath)
  assertthat::assert_that(file.exists(inf.facspath))
  dat.facs <- fread(inf.facspath)
  
  cnames.keep <- !colnames(dat.facs) %in% cnames.remove
  dat.facs.mat <- as.data.frame(dat.facs[, ..cnames.keep])
  cellnumbr <- sapply(dat.facs$Well, function(x) paste("cell", WellToCellNumber(x, totalcols = 24, shift = -1), sep = ""))
  
  jrows <- sapply(dat.facs$Well, function(x) WellToRow)
  jcols <- sapply(dat.facs$Well, function(x) WellToColumn)
  
  rownames(dat.facs.mat) <- cellnumbr
  feats.keep <- which(apply(dat.facs.mat, 2, var) > 0)
  dat.facs.mat.filt <- scale(dat.facs.mat[, feats.keep], center = TRUE, scale = TRUE)
  
  dat.facs.umap <- DoUmapAndLouvain(topics.mat = dat.facs.mat.filt, jsettings = jsettings)
  dat.facs.umap$cellplate <- paste(jplate, dat.facs.umap$cell, sep = "_") 
  dat.facs.umap.merge <- left_join(dat.facs.umap, annots.dat, by = "cellplate") %>%
    mutate(fname = inf.facsname)
  
  return(list(dat.facs.umap.merge = dat.facs.umap.merge, dat.facs.mat.filt = dat.facs.mat.filt))
}


jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 30
jsettings[["min_dist"]] <- 0.1
jsettings[["random_state"]] <- 123

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

# remove weird colnames
cnames.remove <- c("V1", "Well", "Tray X (kw)", "Tray Y (kw)", "TIME", "Tray X", "Tray Y", "Trigger Pulse Width", 
                   "Time 1", "Time 2", "Time 3", "Drop Phase", "ROI Bits 1-16", "ROI Bits 17-32", "Classifier Bits", "Sort Result Bits",
                   "Sort Enable Bits")


# Load FACS ---------------------------------------------------------------

ivec <- seq(9)
names(ivec) <- paste("p", ivec, sep = "")

# jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks



jbatch <- "SL1"
indir <- "/Users/yeung/Dropbox/scCHiC_figs/FACS data/FACS data SL sorts/indexfiles"

# infs.lst <- lapply(jmarks, function(jmark){
#   GetFiles(indir, mark = jmark, batch = jbatch)
# }) %>%
#   unlist()

infs.lst <- GetFiles2(indir, batch = jbatch) %>%
  unlist()

assertthat::assert_that(length(infs.lst) > 0)


# fullname <- "20211202Mous BM SL3 04-01 Plate_001_index.csv"
# fullname <- "20211202Mous BM SL3 04-01 09-03 Plate_010_index.csv"
# stringi::stri_extract(fullname, regex = "\\d+-\\d+")
# # library(stringi)


cnames.remove <- c("V1", "Well", "Tray X (kw)", "Tray Y (kw)", "TIME", "Tray X", "Tray Y", "Trigger Pulse Width", 
                   "Time 1", "Time 2", "Time 3", "Drop Phase", "ROI Bits 1-16", "ROI Bits 17-32", "Classifier Bits", "Sort Result Bits",
                   "Sort Enable Bits")


dat.facs.mat.filt.lst <- lapply(infs.lst, function(jinf){
  dat.facs <- data.table::fread(file = jinf)
  annots.lst <- GetAnnotFromName(jinf)
  cnames.keep <- !colnames(dat.facs) %in% cnames.remove
  dat.facs.mat <- as.data.frame(dat.facs[, ..cnames.keep])
  cellnumbr <- sapply(dat.facs$Well, function(x) paste("cell", WellToCellNumber(x, totalcols = 24, shift = 0), sep = ""))
  wells <- dat.facs$Well
  
  jrows <- sapply(dat.facs$Well, WellToRow)
  jcols <- sapply(dat.facs$Well, WellToColumn)
  
  rownames(dat.facs.mat) <- cellnumbr
  feats.keep <- which(apply(dat.facs.mat, 2, var) > 0)
  dat.facs.mat.filt <- as.data.frame(scale(dat.facs.mat[, feats.keep], center = TRUE, scale = TRUE))
  # dat.facs.mat.filt$Batch <- annots.lst$batch
  # dat.facs.mat.filt$Mark <- annots.lst$mark
  # dat.facs.mat.filt$Plate <- annots.lst$plate
  dat.annot <- data.frame(cellnumbr = cellnumbr, Well = wells, Row = jrows, Column = jcols, Batch = annots.lst$batch, Mark = annots.lst$mark, Plate = annots.lst$plate, stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(cell = paste(cellnumbr, Batch, Mark, Plate, sep = "_"))
  rownames(dat.facs.mat.filt) <- dat.annot$cell
  return(list(mat = dat.facs.mat.filt, meta = dat.annot))
})  
cnames.common <- Reduce(f = intersect, lapply(dat.facs.mat.filt.lst, function(x) colnames(x$mat)))

dat.facs.mat.filt.long <- lapply(dat.facs.mat.filt.lst, function(x){
  jdat <- as.data.frame(x$mat[, cnames.common])
}) %>%
  bind_rows()

dat.meta <- lapply(dat.facs.mat.filt.lst, function(x){
  jdat <- as.data.frame(x$meta)
}) %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(ctype = BatchColumn2Ctype(Batch, Column, Row))
  

rownames(dat.facs.mat.filt.long) <- dat.meta$cell

# do PCA first
pca.out <- prcomp(t(dat.facs.mat.filt.long), center = TRUE, scale. = TRUE)
dat.pca <- pca.out$rotation

# dat.facs.umap <- DoUmapAndLouvain(topics.mat = dat.facs.mat.filt.long, jsettings = jsettings)
dat.facs.umap <- DoUmapAndLouvain(topics.mat = dat.pca, jsettings = jsettings) 

ggplot(dat.facs.umap, aes(x = umap1, y = umap2, color = louvain)) + 
  geom_point() + 
  ggtitle(jbatch) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

BatchColumn2Ctype <- function(Batch, Column, Row = NA){
  
  if (Batch == "SL1"){
    if (Column >= 1 & Column <= 4){
      ctype <- "TCells"
    } else if (Column >= 5 & Column <= 8){
      ctype <- "BCells"
    } else if (Column >= 9 & Column <= 12){
      ctype <- "NKCells"
    } else if (Column >= 13 & Column <= 16){
      ctype <- "Eryths"
    } else if (Column >= 17 & Column <= 24)
      ctype <- "AllCells"
  }
  
  if (Batch == "SL2"){
    if (Column >= 1 & Column <= 4){
      ctype <- "Neutrophils"
    } else if (Column >= 5 & Column <= 8){
      ctype <- "Monocytes"
    } else if (Column >= 9 & Column <= 12){
      ctype <- "DCs"
    } else if (Column >= 13 & Column <= 24){
      ctype <- "All"
    }
  }
  
  if (Batch == "SL3"){
    if (Column == 1){
      ctype <- "HSCs"
    } else if (Column == 2){
      ctype <- "MPPs"
    } else if (Column == 3){
      ctype <- "LT"
    } else if (Column == 4){
      ctype <- "ST"
    } else if (Column >= 5 & Column <= 24){
      ctype <- "LSK"
    } else {
      warning("Error: unknown column:", Column)
    }
  }
  
  if (Batch == "SL4"){
    if (Column == 1){
      ctype <- "GMP"
    } else if (Column == 2){
      ctype <- "CMP"
    } else if (Column == 3){
      ctype <- "MEP"
    } else if (Column >= 4 & Column <= 24){
      ctype <- "LSK"
    }
  }
    
  if (Batch == "SL5"){
    if (Column == 1 & Row >= 1 & Row <= 8){
      ctype <- "pDCs"
    } else if (Column == 1 & Row >= 9 & Row <= 16){
      ctype <- "IL7RLinNeg"
    } else if (Column >= 2 & Column <= 12){
      ctype <- "LinNeg"
    } else if (Column >= 13 & Column <= 24){
      ctype <- "LSK"
    }
  }
    
  return(ctype)
}

dat.facs.umap.annot <- dat.facs.umap %>%
  left_join(., dat.meta) 
  

jmarks.sl <- unique(dat.facs.umap.annot$Mark)

for (jmark in jmarks.sl){
  m <- ggplot(dat.facs.umap.annot %>% filter(Mark == jmark), aes(x = umap1, y = umap2, color = ctype)) + 
    geom_point() + 
    ggtitle(paste(jbatch, jmark)) + 
    facet_grid(Plate~ctype) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}


# Get column and row numbers ----------------------------------------------





