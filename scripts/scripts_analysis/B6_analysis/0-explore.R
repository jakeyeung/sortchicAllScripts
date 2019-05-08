# Jake Yeung
# Date of Creation: 2019-05-06
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/1-remove_bad_bins_define_good_cells.R
# Filter out bad cells and bins 

library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
library(Matrix)
library(hash)
library(parallel)

GetRepB6 <- function(x){
  return(strsplit(x, "_")[[1]][[5]])
}

GetCellB6 <- function(x){
  return(strsplit(x, "_")[[1]][[6]])
}

ParseFilenameB6 <- function(fname){
  # make into H3K4me1.BM.m1.S9.AH3VGVBGX9.CAGAATAT so we can use MARA loading scripts to convert names
  # fname format: PZ-BM-m2-H3K27me3-2_H2GV2BGX9_S18_L002_R1_001.summary.out
  fbase <- strsplit(fname, split = "\\.")[[1]][[1]]
  mark <- strsplit(fbase, split = "-")[[1]][[4]]
  mouse <- strsplit(fbase, split = "-")[[1]][[3]]
  tiss <- strsplit(fbase, split = "-")[[1]][[2]]
  techrep <- strsplit(fbase, split = "_")[[1]][[3]]
  techrepname <- strsplit(fbase, split = "_")[[1]][[2]]
  cellname <- paste(mark, tiss, mouse, techrep, techrepname, sep = ".")  # no cell barcode need to add it yourself 
  return(cellname)
}

ReadAndParseFileB6 <- function(inf, cellhash.bc){
  dat <- data.table::fread(inf, header = FALSE, stringsAsFactors = FALSE, col.names = c("cellname", "bc", "dinuc"))
  dat$cellname.bc <- sapply(dat$bc, function(x) ifelse(!is.null(cellhash.bc[[x]]), cellhash.bc[[x]], NA))
  dat.sum <- subset(dat, !is.na(cellname.bc))%>%
    group_by(cellname.bc, dinuc) %>%
    summarise(dinuc.freq = length(dinuc)) %>%
    group_by(cellname.bc) %>%
    mutate(dinuc.frac = dinuc.freq / sum(dinuc.freq))
  return(dat.sum)
}

FilterAndSwitchColnamesB6 <- function(dat.sum){
  bcs.uniq <- unique(dat.sum$bc)
  bcs.keep <- sapply(bcs.uniq, function(x){
    ifelse(!is.null(cellhash[[x]]), cellhash[[x]], NA)
  })
  bcs.keep <- names(bcs.keep[!is.na(bcs.keep)])
  dat.sum <- subset(dat.sum, bc %in% bcs.keep)
  dat.sum$cname.new <- SwitchColnames(dat.sum$cellname)
  return(dat.sum)
}

# Load TA frequencies -----------------------------------------------------

infs <- list.files("/Users/yeung/data/scchic/from_cluster/count_summaries_B6/parsed2_filt", pattern = "*.out.gz", full.names = TRUE)
barcodes <- read.table("data/barcode_summaries/barcodes/maya_384NLA.bc", col.names = "Barcodes", stringsAsFactors = FALSE)
cellhash <- hash(rownames(barcodes), unlist(barcodes))
cellhash.bc <- hash(unlist(barcodes), paste("cell", rownames(barcodes), sep = ""))

infs <- infs[[1]]
system.time(
  dat.sum <- mclapply(infs, function(inf){
    return(ReadAndParseFileB6(inf, cellhash.bc = cellhash.bc))
  }, mc.cores = 1) %>% 
    bind_rows()
)

ggplot(dat.sum %>% filter(dinuc == "TA"), aes(x = cellname.bc, y = dinuc.frac)) + geom_point() 

# Load data ---------------------------------------------------------------

indir <- "/Users/yeung/data/scchic/from_cluster/count_mat_B6_merged"

jmark <- "H3K4me1"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
names(jmarks) <- jmarks

dat.mat.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- file.path(indir, paste0("BM-B6-merged-", jmark, ".txt"))
  assertthat::assert_that(file.exists(inf))
  dat.mat <- as.data.frame(fread(inf))
  rownames(dat.mat) <- dat.mat$V1
  dat.mat$V1 <- NULL
  dat.mat <- Matrix(as.matrix(dat.mat), sparse=TRUE)
  return(dat.mat)
}) 


# Cell sums ---------------------------------------------------------------


dat.cellsizes <- lapply(jmarks, function(jmark){
  dat.mat <- dat.mat.lst[[jmark]]
  cellsizes <- Matrix::colSums(dat.mat)
  jout <- data.frame(cellsize = cellsizes, mark = jmark, bin = names(cellsizes))
  return(jout)
}) %>%
  bind_rows()

dat.cellsizes$cellname <- sapply(dat.cellsizes$bin, GetCellB6)

indx.all <- seq(384)
hascell.indx <- c(seq(1:356),seq(360:379)+360)
empty.indx <- setdiff(indx.all, hascell.indx)

empty.names <- paste("cell", empty.indx, sep = "")

dat.cellsizes <- dat.cellsizes %>%
  rowwise() %>%
  mutate(is.empty = cellname %in% empty.names)

ggplot(dat.cellsizes, aes(x = cellsize, fill = is.empty)) + geom_density(alpha = 0.25) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark) + scale_x_log10()
ggplot(dat.cellsizes, aes(x = cellsize, fill = is.empty)) + geom_density(alpha = 0.25) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark) + scale_x_log10()

dat.cellsizes.sum <- dat.cellsizes %>%
  group_by(mark) %>%
  summarise(xintercept = median(cellsize))

ggplot(dat.cellsizes, aes(x = cellsize)) + geom_density() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark, ncol = 1)  + geom_vline(data = dat.cellsizes.sum, mapping = aes(xintercept = xintercept)) + scale_x_log10() 


# Take mean across cells fore ach bin -------------------------------------

dat.mean.long <- lapply(jmarks, function(jmark){
  dat.mat <- dat.mat.lst[[jmark]]
  jmeans <- Matrix::rowMeans(dat.mat)
  # make data frame
  jout <- data.frame(cellmean = jmeans, mark = jmark, bin = names(jmeans))
}) %>%
  bind_rows() %>%
  group_by(mark) %>%
  mutate(binrank = rank(cellmean, ties.method = "max"), 
         binrank.norm = binrank / length(binrank),
         cellmean.log = log2(cellmean * 10^6 + 1)) %>%
  arrange(binrank, mark)



# do pairs plot on the cellmean?

dat.wide <- spread(dat.mean.long %>% dplyr::select(mark, bin, cellmean.log), key = "mark", value = "cellmean.log")
dat.rank.wide <- spread(dat.mean.long %>% dplyr::select(mark, bin, binrank), key = "mark", value = "binrank")
dat.rank.norm.wide <- spread(dat.mean.long %>% dplyr::select(mark, bin, binrank.norm), key = "mark", value = "binrank.norm")

library(GGally)

ggpairs(dat.wide %>% dplyr::select(-bin), lower = list(continuous = wrap("points", alpha = 0.2))) + theme_classic()
# ggpairs(dat.rank.norm.wide %>% dplyr::select(-bin), lower = list(continuous = wrap("points", alpha = 0.2))) + theme_classic()

# clip top 1 percent, bottom 1 percent, replot
rank.mean <- data.frame(bin = dat.rank.norm.wide$bin, mark.mean = rowMeans(dat.rank.norm.wide %>% dplyr::select(-bin)))

# bad.bins <- subset(rank.mean, mark.mean >= 0.99 | mark.mean <= 0.01)$
# 
# # Bins at top are bad, bins at bottom also bad ----------------------------
# 
# 
# 
# 
# # Check colnames ----------------------------------------------------------
# 
# # cellnames <- colnames(dat.mat)[2:ncol(dat.mat)]
# 
# # reps <- sapply(cellnames, GetRepB6, simplify = TRUE, USE.NAMES = FALSE)
# # cells <- sapply(cellnames, GetCellB6, simplify = TRUE, USE.NAMES = FALSE)


