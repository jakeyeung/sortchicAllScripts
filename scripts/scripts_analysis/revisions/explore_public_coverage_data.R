# Jake Yeung
# Date of Creation: 2019-08-12
# File: ~/projects/scchic/scripts/scripts_analysis/revisions/explore_public_coverage_data.R
# Explore coverage data

rm(list=ls())

library(JFuncs)
library(dplyr)
library(ggplot2)
library(data.table)
library(topicmodels)
library(tidytext)
library(umap)
library(ggrepel)

library(tidyr)

library(hash)
library(igraph)

library(ChIPseeker)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

library(DropletUtils)
# library(here)



source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/PlotFunctions.R")


# Load data ---------------------------------------------------------------

indir <- "/Users/yeung/data/scchic/from_cluster/revision_objs/public_data_objs"

infs <- list.files(indir, pattern = "*.tab", full.names = TRUE)
ctypes <- sapply(infs, function(inf) strsplit(basename(inf), split = "_")[[1]][[2]], USE.NAMES = FALSE)
names(infs) <- ctypes

dat.chipseq <- lapply(ctypes, function(ctype){
  dat.tmp <- fread(infs[[ctype]])
  dat.tmp$ctype <- ctype
  # take first 4 columns for now, rest are chic
  dat.tmp <- dat.tmp[, c(1:4)]
  colnames(dat.tmp) <- c("chr", "start", "end", ctype)
  return(dat.tmp)
}) %>%
  purrr::reduce(left_join)
# https://stackoverflow.com/questions/8091303/simultaneously-merge-multiple-data-frames-in-a-list/34393416#34393416


# What is the correlation between chip-seq?  ------------------------------

# remove lowly expressed bins and do PCA
jcutoff <- 5

cols.keep <- !colnames(dat.chipseq) %in% c("chr", "start", "end")
dat.filt <- dat.chipseq[which(rowMeans(dat.chipseq[, cols.keep]) > jcutoff), ]

print(dim(dat.filt))
# dat.filt <- dat.chipseq

# take variable bins?


dat.input <- as.matrix(t(subset(dat.filt, select = c(-chr, -start, -end))))

# normalize matrix??
library(scchicFuncs)
library(ggrepel)

# dat.input.norm <- t(NormalizeMatrix(t(dat.input)))
dat.input.norm <- t(scale(t(dat.input), center = TRUE, scale = TRUE))
dat.input.norm <- scale(dat.input.norm, center = TRUE, scale = TRUE)

dat.pca <- prcomp(dat.input.norm, center = FALSE, scale. = FALSE)

dat.proj <- dat.input %*% dat.pca$rotation %*% diag(dat.pca$sdev)

dat.proj.long <- data.frame(celltype = rownames(dat.proj), pc1 = dat.proj[, 1], pc2 = dat.proj[, 2], stringsAsFactors = FALSE)

ggplot(dat.proj.long, aes(x = pc1, y = pc2, label = celltype)) + geom_point() + geom_text_repel() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

plot(dat.pca$x[, 1], dat.pca$x[, 2], pch = 20)
text(dat.pca$x[, 1], dat.pca$x[, 2], labels = rownames(dat.pca$x))


# Normalize magic ---------------------------------------------------------

# handle the NaNs
dat.mat <- as.matrix(dat.chipseq[, -c(1,2,3)])
rownames(dat.mat) <- paste(dat.chipseq$chr, paste(dat.chipseq$start, dat.chipseq$end, sep = "-"), sep = ":")

# set NaNs to zeros (rowwise??)
dat.mat[which(is.nan(dat.mat))] <- 0

# remove the top outlier in each
for (col.i in seq(ncol(dat.mat))){
  print(col.i)
  vec.tmp <- dat.mat[, col.i]
  row.i <- which(vec.tmp == max(vec.tmp))
  dat.mat[row.i, col.i] <- 0
}

boxplot(dat.mat)

plot(density(unlist(dat.mat)))

cnames.tmp <- colnames(dat.mat)
rnames.tmp <- rownames(dat.mat)
dat.mat <- preprocessCore::normalize.quantiles(dat.mat, copy = TRUE)
colnames(dat.mat) <- cnames.tmp
rownames(dat.mat) <- rnames.tmp


boxplot(dat.mat)

# Do multinomial likelihoods ----------------------------------------------

# set up the proportions.... maybe just fractional of log? 





# dat.filt <- dat.mat[apply(dat.mat, 1, function(jrow) !any(is.nan(jrow))), ]

# undo the log transform? 
# dat.mat <- 2^dat.mat

# dat.mat.filt <- sweep(dat.mat, MARGIN = 2, STATS = colSums(dat.mat), FUN = "/")
# dat.mat.filt <- sweep(dat.mat, MARGIN = 2, STATS = colSums(dat.mat), FUN = "/")
dat.mat.filt <- dat.mat

# first do B, EryA, GN, NK, and LT
ctypes.filt <- c("B", "EryB", "GN", "NK", "LT")

probs.lst <- as.list(as.data.frame(dat.mat.filt))
# name the list just to be safe
probs.lst <- lapply(probs.lst, function(x){
  names(x) <- rownames(dat.mat.filt)
  return(x)
}) 


# Calculate lkelihoods ----------------------------------------------------

# get the raw data

# load trajs
inf.traj.stringent <- "/Users/yeung/data/scchic/robjs/B6_objs_stringent/traj_objs_H3K4me3.stringent.Rdata"
load(inf.traj.stringent, v=T)
dat.umap.long.trajs.stringent <- dat.umap.long.trajs$H3K4me3
trajs.stringent <- trajs$H3K4me3
trajs.objs.stringent <- trajs.objs$H3K4me3
inf.traj <- "/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata"
load(inf.traj, v=T)
dat.umap.long.trajs$H3K4me3 <- dat.umap.long.trajs.stringent
trajs$H3K4me3 <- trajs.stringent
trajs.objs$H3K4me3 <- trajs.objs.stringent


jmark <- "H3K4me1"
inf.raw <- paste0("/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA/B6_", jmark, "_pcutoff_0.95_binfilt_cellfilt.2019-05-11.RData")
load(inf.raw, v=T)

bins.keep <- rownames(dat.mat.filt)
# remove chr
# bins.keep <- sapply(bins.keep, function(x) gsub("^chr", "", x), USE.NAMES = FALSE)

# filter for bins 
# add chr to rownames
rownames(count.dat$counts) <- sapply(rownames(count.dat$counts), function(x) paste0("chr", x), USE.NAMES = FALSE)
count.filt <- count.dat$counts[bins.keep, ]

# binarize??
count.filt.bin <- BinarizeMatrix(count.filt)

# check we didnt' filter out too many
par(mfrow=c(2,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
plot(hist(Matrix::colSums(count.filt), breaks = 100))
plot(hist(Matrix::colSums(count.dat$counts), breaks = 100))
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)

# get likelihood

# get a cell from granulocyte, B-cell, and erythro for H3K4me1


# get trajs
cell.i <- 1
cell.gn <- (trajs[[jmark]][["granu"]] %>% arrange(desc(lambda)))$cell[[cell.i]]
cell.b <- (trajs[[jmark]][["lymphoid"]] %>% arrange(desc(lambda)))$cell[[cell.i]]
cell.eryth <- (trajs[[jmark]][["eryth"]] %>% arrange(desc(lambda)))$cell[[cell.i]]
cell.nk <- (trajs[[jmark]][["nk"]] %>% arrange(desc(lambda)))$cell[[cell.i]]
cell.mega <- (trajs[[jmark]][["mega"]] %>% arrange(desc(lambda)))$cell[[cell.i]]
cell.prog <- (trajs[[jmark]][["granu"]] %>% arrange(lambda))$cell[[cell.i]]

# do GN first
LL.lst <- lapply(list(cell.gn, cell.b, cell.eryth, cell.nk, cell.mega, cell.prog), function(cell.name){
  cell.vec <- count.filt[, cell.name]
  LL.vec <- sapply(probs.lst, function(jprob){
    assertthat::assert_that(all(names(cell.vec) == names(jprob)))
    return(dmultinom(x = cell.vec, prob = jprob, log = TRUE))
  })
})

p.lst <- lapply(LL.lst, function(LL){
  SoftMax(LL)
})

print(LL.lst)
print(p.lst)

# Why is EryA not behaving? ------------------------------------------------

# couny.filt.norm <- NormalizeMatrix(count.filt, use.sweep = TRUE)

eryth.vec <- count.filt[, cell.eryth]

# find sox6
eryth.vec.chr7 <- eryth.vec[grepl("^chr7:115", names(eryth.vec))]

# where is sox6 in the proportions vector?
eryth.pvec <- probs.lst$EryA
eryth.pvec.chr7 <- eryth.pvec[grepl("^chr7:115", names(eryth.pvec))]


# Select features for Eryth and redo the LL -------------------------------


jmark <- "H3K4me1"
jbin <- TRUE; kstr <- "25_30_40_50"
inf.lda <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs_bins_B6/ldaAnalysisBins_BinCellFilt2/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.", jbin, ".no_filt/lda_out_meanfilt.B6_", jmark, "_pcutoff_0.CountThres0.K-", kstr, ".Robj")
assertthat::assert_that(file.exists(inf.lda))
kchoose <- "50"
out.objs <- LoadLDABins(jmark, jbin = NA, top.thres = 0.995, inf = inf.lda, convert.chr20.21.to.X.Y = FALSE, add.chr.prefix = TRUE, choose.k = kchoose)


# get probability
# # load data from LDA
indir <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
infs <- lapply(jmarks, function(jmark) list.files(indir, pattern = paste0(jmark, ".RData"), full.names = TRUE))
infs.stringent <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.H3K4me3.stringent_filter.RData"
infs$H3K4me3 <- infs.stringent
tm.result.lst <- lapply(infs, LoadGetTmResult)

# get top hits for the eryth topic
topic.mat <- tm.result.lst[[jmark]]$topics
term.mat <- tm.result.lst[[jmark]]$terms
# gues toipc from sorting by topic weight
# eryth.topic <- names(sort(topic.mat[cell.eryth, ], decreasing=TRUE))[[1]]

# which cell has highest eryth topic?


eryth.topic <- 7
lymph.topic <- 14
lymph.topic2 <- 28
granu.topic <- 47
granu.topic2 <- 3
nk.topic <- 40
mf.topic <- 13

cell.eryth2 <- names(head(sort(topic.mat[, eryth.topic], decreasing = TRUE), n = 50))[[1]]

topics.keep <- c(7,14,28,47,3,40,13)
names(topics.keep) <- topics.keep

nterms.keep <- 1000
terms.keep <- lapply(topics.keep, FUN = function(jtop){
  eryth.terms <- term.mat[jtop, ]
  jterms <- names(head(sort(eryth.terms, decreasing=TRUE), n = nterms.keep))
  return(jterms)
}) %>%
  unlist()

# hande duplicates
terms.keep <- unique(terms.keep)

dat.topics.filt <- dat.mat[unlist(terms.keep), ]
dat.topics.filt <- sweep(dat.topics.filt, MARGIN = 2, STATS = colSums(dat.topics.filt), FUN = "/")

# handle zeros
zero.fill <- min(dat.topics.filt[which(dat.topics.filt > 0)])
dat.topics.filt[which(dat.topics.filt == 0)] <- zero.fill

# plot(density(unlist(dat.topics.filt)))
# abline(v = 1.18e-5)

# check EryA and B top hits

top.ery.hits <- rownames(head(dat.topics.filt[order(dat.topics.filt[, "EryA"], decreasing = TRUE), ], n = 50))

# chekc weight son eryth topic
jsub <- term.mat[, top.ery.hits]
barplot(rowSums(jsub))  # enrichment in 7 found

# let's do the LL now

probs.lst.filt <- as.list(as.data.frame(dat.topics.filt))
# name the list just to be safe
probs.lst.filt <- lapply(probs.lst.filt, function(x){
  names(x) <- rownames(dat.topics.filt)
  return(x)
}) 

# add a fake one to check eryth works
eryth.weights.real <- sort(term.mat[eryth.topic, ])
# rearrange to match rownames
eryth.weights.real <- eryth.weights.real[rownames(dat.topics.filt)]

# probs.lst.filt$eryth.from.LDA <- eryth.weights.real

LL.filt.lst <- lapply(list(cell.gn, cell.b, cell.eryth, cell.nk, cell.mega, cell.prog, cell.eryth2), function(cell.name){
  cell.vec <- count.filt[terms.keep, cell.name]
  print(sum(cell.vec))
  LL.vec <- sapply(probs.lst.filt, function(jprob){
    assertthat::assert_that(all(names(cell.vec) == names(jprob)))
    return(dmultinom(x = cell.vec, prob = jprob, log = TRUE))
  })
})

p.filt.lst <- lapply(LL.filt.lst, SoftMax)

# check the top proportions in EryA and see if it matches the raw counts

eryA.terms <- names(head(sort(probs.lst.filt$EryA, decreasing = TRUE), n = 100))
eryth.terms <- names(head(sort(term.mat[eryth.topic, ], decreasing = TRUE), n = 100))

cell.vec <- count.filt[terms.keep, cell.eryth]

cell.vec[eryA.terms]
cell.vec[eryth.terms]

# Do for all cells --------------------------------------------------------
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)

# sort by louvains maybe

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0", "#32CD32", "#D3D3D3")
PlotXYWithColor(dat.umap.long.trajs[[jmark]], xvar = "umap1", yvar = "umap2", cname = "louvain", cont.color = FALSE, col.palette = cbPalette, jsize = 4)

# do all cells in a louvain
jlouv <- 2  # eryth
jlouv <- 11  # Bcell
# eryth.cells <- subset(dat.umap.long.trajs[[jmark]], louvain == jlouv)$cell
eryth.cells <- subset(dat.umap.long.trajs[[jmark]])$cell
names(eryth.cells) <- eryth.cells

counts.final <- 10000
# counts.final <- 5000
cell.counts <- Matrix::colSums(count.filt)  # approximately account for sliding window
prop.vec <- counts.final / cell.counts
# handle >1
prop.vec[which(prop.vec > 1)] <- 1

count.filt.downsamp <- downsampleMatrix(count.filt, prop.vec)

# subsample: 10000 max
LL.ctype.lst <- lapply(eryth.cells, function(cell.name){
  cell.vec <- count.filt[terms.keep, cell.name]
  # subsample?
  #print(sum(cell.vec))
  LL.vec <- sapply(probs.lst.filt, function(jprob){
    assertthat::assert_that(all(names(cell.vec) == names(jprob)))
    return(dmultinom(x = cell.vec, prob = jprob, log = TRUE))
  })
})
p.ctype.lst <- lapply(LL.ctype.lst, SoftMax)

# summarize into a matrix
jmat <- as.matrix(do.call(rbind, p.ctype.lst))
jmat[is.nan(jmat)] <- -Inf

jdat <- data.frame(cell = rownames(jmat), jmat, stringsAsFactors = FALSE) %>%
  tidyr::gather(key = "celltype", value = "logP", -cell)

ggplot(jdat, aes(x = logP, fill = as.factor(celltype))) + geom_histogram(bins = 50) + facet_wrap(~celltype) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_viridis_d()


# Check why ST shows up  --------------------------------------------------

# downstream, for each cell, pick the best 
cell.names <- names(p.ctype.lst)

cell.counts.downsamp <- Matrix::colSums(count.filt.downsamp)

LL.dat <- lapply(cell.names, function(cname){
  LL.vec <- LL.ctype.lst[[cname]]
  p.vec <- p.ctype.lst[[cname]]
  cell.count = cell.counts[[cname]]
  cell.count.downsamp = cell.counts.downsamp[[cname]]
  if (all(is.infinite(LL.vec))){
    LL.max <- NA
    p.max <- NA
    best.ctype <- NA
  } else {
    LL.max <- max(LL.vec)
    p.max <- max(p.vec)
    best.ctype <- names(which.max(LL.vec))
  }
  dat.tmp <- data.frame(cell = cname, LL.max = LL.max, p.max = p.max, ctype.pred = best.ctype, cell.size = cell.count, cell.count.downsamp = cell.count.downsamp)
  return(dat.tmp) 
}) %>%
  bind_rows()

LL.sum <- LL.dat %>%
  group_by(ctype.pred) %>%
  summarise(ncell = length(ctype.pred))

print(LL.sum)

# plot UMAP with predicted celltype 
LL.dat.merge <- left_join(dat.umap.long.trajs[[jmark]], LL.dat)


cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(LL.dat.merge, 
# ggplot(LL.dat.merge %>% mutate(ctype.pred = ifelse(is.na(ctype.pred), NA, ctype.pred)), 
       aes(x = umap1, y = umap2, color = ctype.pred)) + geom_point(alpha=0.75, size = 3) + 
  scale_color_manual(values = cbPalette) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# collapse EryA and EryB together
# collapse CD4, CD8 together
LL.dat.merge.pretty <- LL.dat.merge %>%
  rowwise() %>%
  mutate(ctype.pred = ifelse(ctype.pred == "EryA" | ctype.pred == "EryB", "EryAorB", ctype.pred),
         ctype.pred = ifelse(ctype.pred == "CD4" | ctype.pred == "CD8", "CD4or8", ctype.pred))

m.pretty <- ggplot(LL.dat.merge.pretty, 
       aes(x = umap1, y = umap2, color = ctype.pred)) + geom_point(alpha=1, size = 3) + 
  scale_color_manual(values = cbPalette) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jmark, "Probabilistic celltype inference using multinomial likelihood, \n with Ido Amit ChIP-seq data as underlying genomic region proportions")


ggplot(LL.dat.merge %>% mutate(ctype.pred = ifelse(is.na(ctype.pred), "znone", ctype.pred)), 
       aes(x = umap1, y = umap2, color = LL.max)) + geom_point() + 
  scale_color_viridis_c() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~ctype.pred)

ggplot(LL.dat %>% mutate(ctype = is.na(ctype.pred)), aes(x = cell.size)) + geom_density() + facet_wrap(~ctype, ncol = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# how diverse are my p-vectors?
plot(density(probs.lst.filt[[1]]))

# Print output
pdf(file = paste0("/Users/yeung/data/scchic/pdfs/compare_with_bulk/multinom_model/", jmark, "_ctype_prediction.", Sys.Date(), ".pdf"), useDingbats = FALSE)
print(m.pretty)
dev.off()


# 
# 
# # Why does this work ? ----------------------------------------------------
# 
# # take a good B cell 
# 
# library(scchicFuncs)
# 
# 
# regions <- data.frame(seqnames = sapply(rownames(count.filt), GetChromo),
#                       start = sapply(rownames(count.filt), GetStart),
#                       end = sapply(rownames(count.filt), GetEnd),
#                       stringsAsFactors = FALSE)
# rownames(regions) <- rownames(count.filt)
# regions <- subset(regions, !seqnames %in% c("chr20", "chr21"))
# 
# regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
# regions.annotated <- as.data.frame(annotatePeak(regions.range,
#                                                 TxDb=TxDb.Mmusculus.UCSC.mm10.knownGene,
#                                                 annoDb='org.Mm.eg.db'))
# regions.annotated$region_coord <- names(regions.range)
# 
# 
# # jcell <- (subset(LL.dat, ctype.pred == ctype) %>% filter(LL.max == min(LL.max)))$cell
# # jcell <- (subset(LL.dat, ctype.pred == ctype) %>% filter(LL.max == min(LL.max)))$cell
# # check erythryoblast
# 
# ctype <- "EryA"
# # ctype <- NA
# # jsub <- subset(LL.dat.merge, umap1 == max(umap1))
# # head(subset(LL.dat.merge, ctype.pred == ctype))
# 
# # debug why znone occurs for so many cells
# jsub <- subset(LL.dat.merge, is.na(ctype.pred)) %>% filter(umap1 == max(umap1))
# jcell <- jsub$cell[[1]]
# 
# 
# # take vector of raw counts
# x <- count.filt[terms.keep, jcell]
# print(sum(x))
# 
# # top hit
# print(head(names(sort(x, decreasing=TRUE))))
# jhit <- names(sort(x, decreasing=TRUE))[[2]]
# 
# jhit <- "chr7:115500000-115600000"  # Sox6
# jhit <- "chr7:115560000-115660000"  # Sox6
# 
# jhit <- "chr7:103720000-103820000"  # Hbb-bt
# subset(regions.annotated, region_coord == jhit)
# 
# 
# # alsod(x) high in the pvec?
# sort(probs.lst.filt[[ctype]], decreasing = TRUE)[1:5]
# 
# probs.lst.filt[[ctype]][names(probs.lst.filt[[ctype]]) == jhit]
# 
# plot(density(probs.lst.filt[[ctype]]))
# abline(v = probs.lst.filt[[ctype]][names(probs.lst.filt[[ctype]]) == jhit])
# 
# # check top hit from reference
# 
# jhits.ref <- names(head(sort(probs.lst.filt[[ctype]], decreasing = TRUE)))
# subset(regions.annotated, region_coord %in% jhits.ref)
# jhit.ref <- names(head(sort(probs.lst.filt[[ctype]], decreasing = TRUE)))[[1]]
# subset(regions.annotated, region_coord == jhit.ref)
# 
# 
# 
# x[grepl(jhit.ref, names(x))]
# 
# 
# plot(density(x))
# abline(v = x[grepl(jhit.ref, names(x))])
# 
# # calculate likelihoods again
# i <- 37
# (jout <- dmultinom(x = x[1:i], prob = probs.lst.filt[[ctype]][1:i], log = TRUE))

