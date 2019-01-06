# Jake Yeung
# Date of Creation: 2019-01-04
# File: ~/projects/scChiC/scripts/scripts_analysis/avo/plot_metacell_knn_and_color.R
# Plot metacell KNN and color

library(hash)
library(ggplot2)
library(dplyr)


# Functions ---------------------------------------------------------------

GetReplHash <- function(experinames.m, biorep, jhash){
  # get replicate names from each mouse
  gstr <- paste0('-BM-', m)
  samps <- grep(gstr, experinames.m, value = TRUE)
  replval.old <- as.numeric(sapply(samps, function(x){
    xtmp <- strsplit(x, "-")[[1]][[4]]
    return(gsub("S", "", xtmp))
  }))
  replval.new <- paste("rep", rank(replval.old), sep = "")
  for (i in seq(length(samps))){
    replhash[samps[i]] <- replval.new[i]
  }
  return(replhash)
}

MakeNewExperiName <- function(old.name, replhash){
  repval <- replhash[[old.name]]
  xshort <- paste(strsplit(old.name, "-")[[1]][c(2,1,3)], collapse = "-")  # switch BM with H3K27me3
  xnew <- paste0(xshort, "-", repval)
  return(xnew)
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

# Cell names --------------------------------------------------------------

barcodes <- read.table("data/barcode_summaries/barcodes/maya_384NLA.bc", col.names = "Barcodes", stringsAsFactors = FALSE)

cellhash <- hash(rownames(barcodes), unlist(barcodes))
cellhash.bc <- hash(unlist(barcodes), paste("cell", rownames(barcodes), sep = ""))


# Plot KNN and note the colors --------------------------------------------



dropboxdir <- "~/Dropbox/scChIC_figs/FIG4_BM"

load(file.path(dropboxdir, "H3K27me3.datadir_mc_f.Rda"))

Q<-object@mc_fp
mc_index<-object@mc
# mc_colors <- sapply(object@mc, function(x) ifelse(x %in% c(5, 7), "#A9A9A9", "#000000"))
mc_colors<-object@colors
colhash <- hash(seq(length(mc_colors)), mc_colors)
cellhash <- hash(names(mc_index), mc_index)
# mc_colors <- rep("#000000", 11)
# mc_colors[c(5, 7)] <- "#A9A9A9"

Q<-data.frame(Q)


load(file.path(dropboxdir, "H3K27me3_2d.datadir_2dproj.Rda"))
x<-object@sc_x
y<-object@sc_y


plot(x,y,pch=16,col=mc_colors[mc_index],cex=.75,axes=FALSE)

# load LDA output 
load("/private/tmp/lda_outputs.meanfilt_1.merge_1000_NoM_binarize.cellmin_1000.cellmax_50000/lda_out.meanfilt.K-20.Robj", v=T)

cellnames <- unname(out.lda[[1]]@documents)
experinames <- unique(sapply(cellnames, function(x) paste(strsplit(x, "-")[[1]][-6], collapse = "-"), USE.NAMES = FALSE))
# add m to replicates
experinames.m <- gsub("-BM-", '-BM-m', experinames)

# rename technical replicates 
replhash <- hash()
mice <- c("m1", "m2")
for (m in mice){
  replhash <- GetReplHash(experinames.m, m, replhash)
}

# add cell barcode
cellnames.new <- sapply(cellnames, MakeNewCellName, replhash, cellhash.bc, add.m=TRUE, USE.NAMES = FALSE)

cellnames.hash <- hash(cellnames.new, cellnames)

# create new mc_index
mc_index.newnames <- mc_index
names(mc_index.newnames) <- sapply(names(mc_index), function(x) cellnames.hash[[x]])



