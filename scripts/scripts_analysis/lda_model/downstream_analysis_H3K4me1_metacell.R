# Jake Yeung
# Date of Creation: 2019-01-22
# File: ~/projects/scChiC/downstream_analysis_H3K4me1_metacell.R
# Recover the Hbb from H3K4me1 analysis from MetaCell

library(topicmodels)
library(dplyr)
library(ggplot2)
library(umap)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(rGREAT)
library(hash)
library(JFuncs)
library(umap)

source("scripts/Rfunctions/MetricsLDA.R")
source("scripts/Rfunctions/AuxLDA.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/GetMetaCellHash.R")


jchip <- "H3K4me1"

jdist <- 1000L
jmean <- 1
jmin <- 100L
jmax <- 500000L
# binarize <- "TRUE";  jtops <- "5_7_10_12_15_20_25_30"
binarize <- "FALSE"; jtops <- "15_20_25_30_35"

jdir <- paste0('/tmp/ldaAnalysisHiddenDomains_', jdist, '/lda_outputs.meanfilt_', jmean, '.cellmin_', jmin, '.cellmax_', jmax, '.binarize.', binarize)
inf <- file.path(jdir, paste0('lda_out_meanfilt.PZ-BM-', jchip, '.CountThres0.K-', jtops, '.Robj'))
infbase <- basename(inf)
infbase <- strsplit(infbase, ".Robj")[[1]][[1]]

# 0.98 threshold 
# inf.GREAT <- file.path(jdir, "downstream", paste0(infbase, ".GREAT.Robj"))
# 0.96 threshold 
inf.GREAT <- file.path(jdir, "downstream", paste0(infbase, ".GREAT.0.96.Robj"))
inf.mc <- file.path(paste0("~/Dropbox/scCHiC_figs/FIG4_BM/", jchip, ".datadir_mc_f.Rda"))
inf.mc.2d <- file.path(paste0("~/Dropbox/scCHiC_figs/FIG4_BM/", jchip, "_2d.datadir_2dproj.Rda"))

inf.mc.mat <- file.path(paste0("/private/tmp/binned_mat/BM-", jchip, "-100kb.txt"))

assertthat::assert_that(file.exists(inf))
assertthat::assert_that(file.exists(inf.GREAT))
assertthat::assert_that(file.exists(inf.mc))
assertthat::assert_that(file.exists(inf.mc.2d))
assertthat::assert_that(file.exists(inf.mc.mat))

mc.out <- GetMCObjects(inf = inf.mc)
load(inf.mc.2d, v=T)
mc_index <- mc.out$mc_index; mc_colors <- mc.out$mc_colors


# plot KNN neighbors
x<-object@sc_x
y<-object@sc_y
plot(x,y,pch=16,col=mc_colors[mc_index],cex=.75,axes=FALSE)

plot(x,y,pch=16,cex=.75,axes=FALSE, col=mc_colors[mc_index])
text(x,y,labels = mc_index)

# check cluster 1

mcell_mc_plot_marks(mc_id=object)

Q <- mc.out$exprs

# diff eprs
allgroup <- sort(unique(mc_index))
ingroup <- c(1)  # 5 is most outer group, 7 is less    # group 5 shows 0.6 fold change. group 7 shows 0.3 FC
outgroup <- allgroup[!allgroup %in% ingroup]

if (length(ingroup) > 1){
  Q<-data.frame(cbind(rowMeans(Q[,c(ingroup)]),rowMeans(Q[,c(outgroup)])))
} else {
  Q<-data.frame(cbind(Q[, ingroup],rowMeans(Q[,c(outgroup)])))
}
Q$fullcoord <- rownames(Q)
Q <- Q %>%
  mutate(chromo = paste("chr", sapply(fullcoord, function(x) strsplit(x, "_")[[1]][[1]]), sep = ""),
         start = sapply(fullcoord, function(x) strsplit(x, "_")[[1]][[2]]),
         end = sapply(fullcoord, function(x) strsplit(x, "_")[[1]][[3]]),
         coord = paste(chromo, paste(start, end, sep = "-"), sep = ":"))

plot(log2(Q[,1]),log2(Q[,2]),pch=16,cex=.25,xlab='average log2-fold change metacell 1',ylab='average log2-fold change other metacells')

# show top hits
head(Q %>% arrange(desc(X1)))


# What is expression of top regions in single cell matrix? ----------------

inmat <- data.table::fread(inf.mc.mat, stringsAsFactors = FALSE)
inmat <- as.data.frame(inmat)
rownames(inmat) <- inmat$V1
inmat$V1 <- NULL

# filter common peaks and cells
peaks.common <- intersect(rownames(inmat), rownames(mc.out$exprs))
cells.common <- intersect(colnames(inmat), names(mc.out$mc_index))

inmat <- inmat[peaks.common, cells.common]
inmat <- sweep(inmat, MARGIN = 2, STATS = colSums(inmat), FUN = "/")

# find this 20% increase in single cells

top.hit <- "chr7:103800000-103900000"
top.hit.format <- gsub("chr", "", top.hit)
top.hit.format <- gsub(":", "_", top.hit.format)
top.hit.format <- gsub("-", "_", top.hit.format)
inmat.sub <- inmat[top.hit.format, ]

dat.out <- data.frame(exprs = unlist(inmat.sub), cell = names(inmat.sub), peak = top.hit)

dat.out$mc <- mc_index[as.character(dat.out$cell)]

ggplot(dat.out, aes(x = as.character(mc), y = log2(exprs + 1))) + 
  geom_jitter(height = 0, width = 0.1) + geom_violin() + 
  stat_summary(aes(y = log2(exprs + 1), group=mc), fun.y=mean, colour="red", geom="line",group=1)
  

