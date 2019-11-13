# Jake Yeung
# Date of Creation: 2019-11-12
# File: ~/projects/scchic/scripts/scripts_analysis/zebrafish/compare_TSS_with_bins_LDA.R
# Compare TSS with bins LDA (variance)


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(hash)
library(igraph)
library(umap)

library(topicmodels)

library(scchicFuncs)

library(ggrepel)

GetGeneAnnotsHash <- function(inf.annot){
  dat.annot <- data.table::fread(inf.annot, col.names = c("chromo", "start", "end", "bname"))
  # add chr
  dat.annot$chromo <- paste("chr", dat.annot$chromo, sep = "")
  rnames.old <- paste(dat.annot$chromo, paste(dat.annot$start, dat.annot$end, sep = "-"), sep = ":")
  rnames.new <- dat.annot$bname
  annots.hash <- hash::hash(rnames.old, rnames.new)
}

AddGeneNameToRows <- function(mat, annots.hash){
  # mat rownmaes got stripped of gene names, add them back
  rnames.old <- rownames(mat)
  rnames.new <- sapply(rnames.old, function(x) annots.hash[[x]])
  rownames(mat) <- rnames.new
  return(mat)
}

# Load data ---------------------------------------------------------------


jmark <- "H3K4me3"
winsize <- 100000L
# winsize <- 50000L
jprefix <- "ZFWKM"
# jprefix <- "ZFWKMCD41plus"
# jprefix <- "ZFWKMCD41plus"

inf.annot <- paste0("/Users/yeung/data/scchic/tables/gene_tss.winsize_", winsize, ".species_drerio.nochr.bed")

# init
inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisTSS_ZFbonemarrow.v2/lda_outputs.PZ-ChIC-", jprefix, "-", jmark, 
              ".winsize_", winsize, ".merged.K-30.binarize.FALSE/ldaOut.PZ-ChIC-", jprefix, "-", jmark, 
              ".winsize_", winsize, ".merged.K-30.Robj")
assertthat::assert_that(file.exists(inf))

x <- load(inf, v=T)

if (length(out.lda) > 1){
  out.lda <- out.lda[[1]]
} 

# add gene name to the coordinates (got lost in mat to sparse mat pipeline)
topics.mat <- posterior(out.lda)$topics
terms.mat <- posterior(out.lda)$terms

colnames(topics.mat) <- paste0("topic_", colnames(topics.mat))
rownames(terms.mat) <- paste0("topic_", rownames(terms.mat))

print(head(out.lda@terms))

# # if no gene names in rownames, then add them 
# annots.hash <- GetGeneAnnotsHash(inf.annot)
# count.mat <- AddGeneNameToRows(count.mat, annots.hash)

# out.lda@terms <- sapply(out.lda@terms, function(x) annots.hash)


# colnames(terms.mat) <- sapply(colnames(terms.mat), function(x) annots.hash[[x]])

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jsettings.1d <- jsettings; jsettings.1d$n_components <- 1

umap.out <- umap(topics.mat, config = jsettings)
umap.out.1d <- umap(topics.mat, config = jsettings.1d)

dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2])
dat.umap.long.1d <- data.frame(cell = rownames(umap.out$layout), umap1.1d = umap.out$layout[, 1])

m.umap.blank <- ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste("LDA with windows around TSS. Winsize:", winsize))

# get variance??
dat.impute.log <- log2(t(topics.mat %*% terms.mat))
rownames(dat.impute.log) <- gsub(";", "_", rownames(dat.impute.log))

# intrachromosomal variance doesnt make sense when doing TSS, do genome wide
# jchromos <- paste("chr", seq(25), sep = "")
jchromos <- c("")

dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.umap.long.varmerge <- left_join(dat.umap.long, dat.var)
dat.umap.long.varmerge <- left_join(dat.umap.long.varmerge, dat.umap.long.1d)

ggplot(dat.umap.long.varmerge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c()

# PlotXYWithColor(dat.umap.long.varmerge, xvar = "umap1", yvar = "umap2", cname = "cell.var.within.sum.norm")

ggplot(dat.umap.long.varmerge, aes(x = cell.var.within.sum.norm, y = umap1.1d)) + geom_point() 



# Compare with bins -------------------------------------------------------


# Load dat ----------------------------------------------------------------

jmark <- "H3K4me3"
jbin <- "FALSE"

# inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZFbonemarrow/lda_outputs.ZF-ZFWKM-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.binarize.", jbin, "/lda_out_meanfilt.ZF-ZFWKM-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.Robj")
# inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZFbonemarrow/lda_outputs.ZF-ZFWKMCD41plus-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.binarize.", jbin, "/lda_out_meanfilt.ZF-ZFWKMCD41plus-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.Robj")
inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZFbonemarrow/lda_outputs.ZF-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.binarize.", jbin, "/lda_out_meanfilt.ZF-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.Robj")
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


# Do variance calculation -------------------------------------------------

jfac <- 10^6
jpseudo <- 0
dat.impute.log <- log2(t(posterior(out.lda[[1]])$topics %*% posterior(out.lda[[1]])$terms) * jfac + jpseudo)
jchromos.num <- seq(25)
jchromos <- paste("chr", jchromos.num, sep = "")

cells.var.chromo.merged <- CalculateVarAll(dat.impute.log, jchromos)


# Plot with variance  -----------------------------------------------------

dat.umap.long.merge <- left_join(dat.umap.long, cells.var.chromo.merged)

ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
  scale_colour_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(jmark)


# Compare bin versus TSS --------------------------------------------------

dat.umap.bin.vs.tss <- left_join(dat.umap.long.merge, dat.umap.long.varmerge, by = "cell")

ggplot(dat.umap.bin.vs.tss, aes(x = umap1.x, y = umap2.x, color = cell.var.within.sum.norm.y)) + geom_point() + theme_bw() + 
  scale_colour_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(jmark)

ggplot(dat.umap.bin.vs.tss, aes(x = umap1.y, y = umap2.y, color = cell.var.within.sum.norm.x)) + geom_point() + theme_bw() + 
  scale_colour_viridis_c() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(jmark)



