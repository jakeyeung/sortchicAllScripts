# Jake Yeung
# Date of Creation: 2020-08-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/explore_index_data.R
# 


rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Functions ---------------------------------------------------------------


DoUmapOnFACS <- function(indir, jplate.i){
  jplate <- paste0("p", jplate.i)
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
  rownames(dat.facs.mat) <- cellnumbr
  feats.keep <- which(apply(dat.facs.mat, 2, var) > 0)
  dat.facs.mat.filt <- scale(dat.facs.mat[, feats.keep], center = TRUE, scale = TRUE)
  
  dat.facs.umap <- DoUmapAndLouvain(topics.mat = dat.facs.mat.filt, jsettings = jsettings)
  dat.facs.umap$cellplate <- paste(jplate, dat.facs.umap$cell, sep = "_")
  dat.facs.umap.merge <- left_join(dat.facs.umap, annots.dat, by = "cellplate") %>%
    mutate(fname = inf.facsname)
  
  return(list(dat.facs.umap.merge = dat.facs.umap.merge, dat.facs.mat.filt = dat.facs.mat.filt))
}





# Load index data ---------------------------------------------------------

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/facs_index_data/K562_cellcycle_2020-08-21/Index"


LoadIndexPlotPCA <- function(inf.check, jtitle, return.unnorm = FALSE){
  dat.check <- fread(inf.check)
  
  # remove weird colnames
  cnames.remove <- c("V1", "Well", "Tray X (kw)", "Tray Y (kw)", "TIME", "Tray X", "Tray Y", "Trigger Pulse Width",
                     "Time 1", "Time 2", "Time 3", "Drop Phase", "ROI Bits 1-16", "ROI Bits 17-32", "Classifier Bits", "Sort Result Bits",
                     "Sort Enable Bits")
  
  cnames.keep <- !colnames(dat.check) %in% cnames.remove
  dat.check.filt <- data.frame(dat.check[, ..cnames.keep], stringsAsFactors = FALSE)
  rownames(dat.check.filt) <- dat.check$Well
  
  # get cellcycle
  colcoords <- as.numeric(sapply(rownames(dat.check.filt), function(x) substr(x = x, start = 2, stop = nchar(x))))
  cc <- sapply(colcoords, function(colcoord) AddCellCycleLabel(colcoord))
  
  cvec <- sapply(cc, function(x){
    if (x == "0_G1"){
      return("red")
    } else if (x == "1_S"){
      return("green")
    } else if (x == "2_G2/M"){
      return("blue")
    }
  })
  # cchash <- hash::hash(rownames())
  if (return.unnorm){
    return(dat.check.filt)
  }
  
  dat.check.filt <- scale(dat.check.filt, center = TRUE, scale = TRUE)
  dat.check.filt <- t(scale(t(dat.check.filt), center = TRUE, scale = TRUE))
  pca.out <- prcomp(dat.check.filt, center = FALSE, scale. = FALSE)
  plot(pca.out$x[, 1], pca.out$x[, 2], col = cvec, main = jtitle)
  return(dat.check.filt)
}

# infs <- list.files(path = indir, pattern = ".*.E+.*.4me3.*.csv", full.names = TRUE)
infs <- list.files(path = indir, pattern = ".*.E+.*.csv", full.names = TRUE)
names(infs) <- sapply(infs, function(x) strsplit(basename(x), "\\.")[[1]][[1]])

print(infs)

dats <- lapply(infs, function(inf){
  bname <- strsplit(basename(inf), "\\.")[[1]][[1]]
  print(bname)
  LoadIndexPlotPCA(inf, bname)
})

dats.unnorm <- lapply(infs, function(inf){
  bname <- strsplit(basename(inf), "\\.")[[1]][[1]]
  print(bname)
  LoadIndexPlotPCA(inf, bname, return.unnorm = TRUE)
})

# check which feature is hoescht staining

inf.test <- infs[["20200703 K562 E+ 27me3_004_index"]]
dat.raw <- fread(inf.test)

dat.test <- dats$`20200703 K562 E+ 27me3_004_index`
dat.test <- scale(dat.test, center = TRUE, scale = TRUE)
dat.test <- t(scale(t(dat.test), center = TRUE, scale = TRUE))

colcoords <- as.numeric(sapply(rownames(dat.test), function(x) substr(x = x, start = 2, stop = nchar(x))))
cc <- sapply(colcoords, function(colcoord) AddCellCycleLabel(colcoord))

cvec <- sapply(cc, function(x){
  if (x == "0_G1"){
    return("red")
  } else if (x == "1_S"){
    return("green")
  } else if (x == "2_G2/M"){
    return("blue")
  }
})

pca.out.test <- prcomp(dat.test, center = TRUE, scale. = TRUE)

plot(pca.out.test$x[, 1], pca.out.test$x[, 2], col = cvec)
text(pca.out.test$rotation[, 1], pca.out.test$rotation[, 2], labels = rownames(pca.out.test$rotation))


sort(pca.out.test$rotation[, 1], decreasing = TRUE)
sort(pca.out.test$rotation[, 1], decreasing = FALSE)

cnames.keep <- c(names(sort(pca.out.test$rotation[, 1], decreasing = TRUE))[[1]], names(sort(pca.out.test$rotation[, 1], decreasing = FALSE))[[1]])


# Explore index data  -----------------------------------------------------

plot(dat.test[, cnames.keep], col = cvec)


cnames.keep <- grep("405", colnames(dat.test), value = TRUE)

plot(dat.test[, cnames.keep[[1]]], dat.test[, cnames.keep[[2]]], col = cvec)


# Find which is hoescht channel  ------------------------------------------

for (cname in cnames.keep){
  plot(density(dat.test[, cname]), main = cname, col = cvec)
}

dat.test.long <- melt(dat.test)
colnames(dat.test.long) <- c("well", "feature", "value")

dat.cvec <- data.table(well = rownames(dat.test), cellcycle = cc, stringsAsFactors = FALSE)

dat.test.long.annot <- dat.test.long %>%
  left_join(., dat.cvec, by = "well")

ggplot(dat.test.long.annot %>% filter(feature %in% cnames.keep), aes(x = value, fill = cellcycle)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~feature, ncol = 4)


# Plot the Hoescht signal  ------------------------------------------------

dat.test <- dats$`20200703 K562 E+ 27me3_004_index`

jnames <- names(dats)
names(jnames) <- jnames

cname <- "X.405..460.50"
dat.hoesch <- lapply(jnames, function(jname){
  jdat <- dats.unnorm[[jname]]
  cname.i <- which(colnames(jdat) == cname)
  print(cname.i)
  jdat.out <- data.frame(well = rownames(jdat), hoesch = jdat[, cname.i], experi = jname, stringsAsFactors = FALSE) %>%
    left_join(., dat.cvec)
}) %>%
  bind_rows()
 
ggplot(dat.hoesch, aes(x = hoesch, fill = cellcycle)) + 
  geom_density(alpha = 0.25)  + 
  facet_wrap(~experi, ncol = 5) + 
  theme_bw(4) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.x = element_text(size = 8))


# Save to output ----------------------------------------------------------

# add col.i and row.i 

WellToRowAndCol <- function(well, rtrn = c("row", "col", "both")){
  # A1 -> row1, column1
  myLetters <- toupper(letters[1:26])
  row.alphabet <- substr(well, start = 1, stop = 1)
  row.i <- match(x = row.alphabet, table = myLetters)
  col.i <- as.numeric(strsplit(well, split = row.alphabet)[[1]][[2]])
  if (rtrn == "both"){
    return(c(row.i, col.i))
  } else if (rtrn == "row"){
    return(row.i)
  } else if (rtrn == "col"){
    return(col.i)
  }
}

dat.hoesch <- dat.hoesch %>%
  rowwise() %>%
  mutate(row.indx = WellToRowAndCol(well, rtrn = "row"),
         col.indx = WellToRowAndCol(well, rtrn = "col"))



# Save to output  ---------------------------------------------------------

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/spikein_fits_cellcycle.hoesch_staining"
outf <- file.path(outdir, paste0("hoesch_on_K562_plates.rds"))
saveRDS(dat.hoesch, file = outf)




