# Jake Yeung
# Date of Creation: 2019-11-13
# File: ~/projects/scchic/scripts/scripts_analysis/zebrafish/zebrafish_FACS_index_analysis.R
# 
rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(topicmodels)

library(hash)
library(igraph)
library(umap)

library(scchicFuncs)

library(xlsx)

library(DESeq2)

library(preprocessCore)

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ChIPseeker)

library(ggrepel)

WellToCellNumber <- function(w, totalcols = 24){
  # Well name (e.g., P20 to Plate coordinate (e.g. 380)
  jlet <- toupper(letters[1:26])
  jrow <- as.numeric(match(gsub("[^a-zA-Z]", "", w), jlet))
  jcol <- as.numeric(gsub("[a-zA-Z]", "", w))
  return( as.character((jrow - 1) * totalcols + jcol ) )
}



# Load dat ----------------------------------------------------------------

jmark <- "H3K9me3"
jbin <- "FALSE"

inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZFbonemarrow/lda_outputs.ZF-ZFWKM-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.binarize.", jbin, "/lda_out_meanfilt.ZF-ZFWKM-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.Robj")
assertthat::assert_that(file.exists(inf))

load(inf, v=T)

# Plot data ---------------------------------------------------------------

kvec <- sapply(out.lda, function(x) x@k)
kchoose <- 30
kchoose.i <- which(kvec == kchoose)
topics.mat <- posterior(out.lda[[kchoose.i]])$topics
colnames(topics.mat) <- paste("topic", colnames(topics.mat), sep = "_")

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(topics.mat, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2]) %>%
  rowwise()
# mutate(is.stem = grepl("CD41plus", cell))

m.umap.first <- ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(jmark)



# Add FACs information  ---------------------------------------------------
# for H3K4me1 only

indir <- "/Users/yeung/data/scchic/facs"

if (jmark == "H3K4me1"){
  gstr <- ".*WKM.*K4.*ME1.*.csv"
} else if (jmark == "H3K4me3"){
  gstr <- ".*WKM.*K4.*ME3.*.csv"
} else if (jmark == "H3K9me3"){
  gstr <- ".*WKM.*K9.*ME3.*.csv"
} else if (jmark == "H3K27me3"){
  gstr <- ".*WKM.*K27Me3*.csv"
}

(flist <- list.files(path = indir, pattern = gstr, full.names = TRUE, recursive = TRUE))



# load files
dat.facs <- fread(flist[[1]])

# remove weird colnames
# cnames.remove <- c("V1", "Well", "Tray X (kw)", "Tray Y (kw)", "TIME")

cnames.remove <- c("V1", "Well", "Tray X (kw)", "Tray Y (kw)", "TIME", "Time 1", "Time 2", "Time 3", "Drop Phase", "Tray X", "Tray Y", "Sort Enable Bits", "Classifier Bits", "ROI Bits 17-32", "ROI Bits 1-16", "Tray X (kw)", "Tray Y (kw)")
cnames.keep <- !colnames(dat.facs) %in% cnames.remove

# cnames.keep <- c("SSC", "FSC")

dat.facs.mat <- as.data.frame(dat.facs[, ..cnames.keep])

rownames(dat.facs.mat) <- dat.facs$Well

# remove bad features
feats.keep <- which(apply(dat.facs.mat, 2, var) > 0)

dat.facs.mat.filt <- scale(dat.facs.mat[, feats.keep], center = TRUE, scale = TRUE)

dat.facs.pca <- prcomp(dat.facs.mat.filt, retx = TRUE, center = FALSE, scale. = FALSE)

plot(dat.facs.pca$x[, 1], dat.facs.pca$x[, 2])

# plot(dat.facs.pca$x[, 1], dat.facs.pca$x[, 2], xlim = c(-10, 10), ylim = c(-10, 10))
# plot(dat.facs.pca$x[, 2], dat.facs.pca$x[, 3])

pca.dist = sqrt(dat.facs.pca$x[, 1] ^ 2 + dat.facs.pca$x[, 2] ^ 2)
pca.dist.dat <- data.frame(well = names(pca.dist), pca.dist = pca.dist, stringsAsFactors = FALSE)

# assign Well to cell number


facs.meta <- data.frame(well = dat.facs$Well, cellindx = sapply(dat.facs$Well, WellToCellNumber), stringsAsFactors = FALSE) %>%
  left_join(., pca.dist.dat)

dat.umap.long.facs <- dat.umap.long %>%
  rowwise() %>%
  mutate(cellindx = strsplit(as.character(cell), "_")[[1]][[2]]) %>%
  left_join(., facs.meta)

# PlotXYWithColor(dat.umap.long.facs %>% dplyr::filter(pca.dist < 20), xvar = "umap1", yvar = "umap2", cname = "pca.dist", cont.color = TRUE) + scale_color_viridis_c()
PlotXYWithColor(dat.umap.long.facs, xvar = "umap1", yvar = "umap2", cname = "pca.dist", cont.color = TRUE) + scale_color_viridis_c()


#

