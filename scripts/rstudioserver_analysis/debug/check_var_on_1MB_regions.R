# Jake Yeung
# Date of Creation: 2020-01-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/debug/check_var_on_1MB_regions.R
# Does 1MB region help? 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)
library(scchicFuncs)

library(parallel)
library(hash)
library(igraph)
library(umap)


# Functions ---------------------------------------------------------------



# xlogp <- Vectorize(FUN = function(xi, pi){
#   if (xi == 0 & pi == 0){
#     return(0)
#   } else {
#     xi * log(pi)
#   }
# }, vectorize.args = c("xi", "pi"))

xlogp <- Vectorize(function(xi, pi){
  if (xi == 0 & pi == 0){
    return(0)
  } else {
    xi * log(pi)
  }
}, vectorize.args = c("xi", "pi"))

multinomdev <- function(x, p){
  -2 * sum(xlogp(x, p))
  # -2 * xlogp(x, p)
}

multinomdevdiff <- function(xvec, pnull = NULL, softmax.out = TRUE, remove.zeros = TRUE){
  if (remove.zeros){
    xvec <- xvec[which(xvec > 0)] 
  }
  # dev.sat <- multinomdev(x = xvec, p = xvec / sum(xvec))
  dev.sat <- multinomial_deviance(x = xvec, xvec / sum(xvec))
  if (is.null(pnull)){
    pnull <- rep(1 / length(xvec), length(xvec))
  }
  # dev.null <- multinomdev(x = xvec, p = pnull)
  dev.null <- multinomial_deviance(x = xvec, p = pnull)
  if (softmax.out){
    SoftMax(c(dev.null - dev.sat, dev.null), return.log = TRUE)
  } else {
    return(dev.null - dev.sat)
  }
}

# Load data ---------------------------------------------------------------

inf <- "/home/jyeung/data/from_rstudioserver/intestinalchic/var_across_intestinal_marks.rds"

out.lst <- readRDS(inf)


inf.test <- "/home/jyeung/hpc/intestinal_scchic/raw_data/HVG-intestines.tagged_bams.all_merged/countTables_otherWinSize/HVG_Scraped_Bl6_k4me3_2019-12-20.1000000.1000000.countTable.csv"
mat <- ReadMatSlideWinFormat(inf.test, as.sparse = TRUE)

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")


# chr.keep <- "^chr1:"
# chr.keep <- "^chr19"
chr.keep <- ""
terms.keep <- grep(pattern = chr.keep, x = rownames(mat), value = TRUE)
# order by start site
terms.keep.sorted <- terms.keep[order(sapply(terms.keep, function(x) as.numeric(GetStart(x))))]
starts <- sapply(terms.keep, function(x) as.numeric(GetStart(x)))
plot(starts[terms.keep.sorted])

cells.keep <- out.lst$k4me3$dat.meta$cell

cells.keep.high <- subset(out.lst$k4me3$dat.meta, ncuts > 10000)$cell
mat.filt <- mat[terms.keep, cells.keep]

dat.ncuts <- data.frame(cell = colnames(mat.filt), ncuts2 = colSums(mat.filt), 
                        sparsity = apply(mat.filt, 2, function(xvec) 1 - nnzero(xvec) / length(xvec)), 
                        stringsAsFactors = FALSE)

ggplot(dat.ncuts, aes(x = ncuts, y = sparsity)) + geom_point()

# plot(mat.filt[, 1])
# plot(mat.filt[, 2])
# plot(mat.filt[, 3])
# plot(mat.filt[, 4])

head(out.lst$k4me3$dat.meta)

# get deviance for each cell

system.time(
  dev.diff.lst <- apply(mat.filt, 2, multinomdevdiff, pnull = rep(1 / nrow(mat.filt), nrow(mat.filt)), softmax.out = FALSE)
)

dev.diff.dat <- data.frame(cell = names(dev.diff.lst), devdiff = dev.diff.lst, stringsAsFactors = FALSE) %>%
  left_join(., out.lst$k4me3$dat.meta) %>%
  left_join(., dat.ncuts)

ggplot(dev.diff.dat, aes(y = devdiff, x = ncuts, color = cell.var.within.sum.norm)) + geom_point() + 
  # scale_x_log10() + scale_y_log10() + 
  scale_color_viridis_c(direction = -1) 

ggplot(dev.diff.dat, aes(x = devdiff, y = cell.var.within.sum.norm, color = log10(ncuts))) + geom_point() + 
  scale_color_viridis_c(direction = 1) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Deviance\n(saturated multinom model vs null multinom model)")

ggplot(dev.diff.dat, aes(x = devdiff, y = cell.var.within.sum.norm, color = log10(ncuts))) + geom_point() + 
  scale_color_viridis_c(direction = 1) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Deviance\n(saturated multinom model vs null multinom model)")

ggplot(dev.diff.dat, aes(x = devdiff, y = log10(ncuts), color = cell.var.within.sum.norm)) + geom_point() + 
  scale_color_viridis_c(direction = -1) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
  xlab("Deviance\n(saturated multinom model vs null multinom model)")

ggplot(dev.diff.dat, aes(x = devdiff, y = ncuts, color = cell.var.within.sum.norm)) + geom_point() + 
  scale_color_viridis_c(direction = -1) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
  xlab("Deviance\n(saturated multinom model vs null multinom model)") + 
  scale_y_log10()

ggplot(dev.diff.dat, aes(y = sparsity, x = ncuts, color = cell.var.within.sum.norm)) + geom_point() + 
  scale_x_log10() + scale_y_log10() +
  scale_color_viridis_c(direction = -1) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# visualize a low and high deviance cell
cvec <- sapply(rownames(mat.filt), GetChromo)
cvec.col <- sapply(cvec, function(x) ifelse(x %in% c("chr11", "chr19"), "blue", "red"))

ncuts.min <- 10000
ncuts.max <- 12000
# high.var.cell <- (out.lst$k4me3$dat.meta %>% filter(ncuts > ncuts.min & ncuts < ncuts.max ) %>% arrange(desc(cell.var.within.sum.norm)))$cell[[1]]
# low.var.cell <- (out.lst$k4me3$dat.meta %>% filter(ncuts > ncuts.min & ncuts < ncuts.max ) %>% arrange(cell.var.within.sum.norm))$cell[[1]]

(high.var.cell.dat <- (subset(dev.diff.dat %>% filter(ncuts > ncuts.min & ncuts < ncuts.max)) %>% arrange(desc(cell.var.within.sum.norm)))[1, ])
(low.var.cell.dat <- (subset(dev.diff.dat %>% filter(ncuts > ncuts.min & ncuts < ncuts.max)) %>% arrange(cell.var.within.sum.norm))[1, ])

high.var.cell <- high.var.cell.dat$cell
low.var.cell <- low.var.cell.dat$cell

xvec.high <- mat.filt[, high.var.cell]
xvec.low <- mat.filt[, low.var.cell]

par(mfrow=c(2,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
  plot(xvec.high, main = paste("High Var Cell:", signif(high.var.cell.dat$devdiff, digits = 2), 
                               signif(high.var.cell.dat$cell.var.within.sum.norm, digits = 2), 
                               signif(high.var.cell.dat$ncuts, digits = 2)), 
       col = cvec.col, type = "o")
  plot(xvec.low, main = paste("Low Var Cell:", signif(low.var.cell.dat$devdiff, 2), 
                              signif(low.var.cell.dat$cell.var.within.sum.norm, 2), 
                              signif(low.var.cell.dat$ncuts, 2)), 
       col = cvec.col, type = "o")
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)


# 
# 
# xvec.high <- mat.filt[, high.var.cell]
# xvec.low <- mat.filt[, low.var.cell] 
# 
# dev.diff.high <- multinomdevdiff(xvec.high)
# dev.diff.low <- multinomdevdiff(xvec.low)
# 
# xvec.tmp <- mat.filt[, "HVG-Black6-intestine-k4me3-index4-06112019_90"]
# xvec.tmp <- xvec.low
# 
# dev.sat <- multinomdev(x = xvec.tmp, p = xvec.tmp / sum(xvec.tmp))
# dev.null <- multinomdev(x = xvec.tmp, p = rep(1 / length(xvec.tmp), length(xvec.tmp)))
# (dev.diff <- dev.null - dev.sat)
# 1 - pchisq(dev.diff, length(xvec.tmp) - 1)
# SoftMax(c(dev.null - dev.sat, dev.null))
# 
# icv.high <- signif(var(xvec.high), digits = 2)
# icv.low <- signif(var(xvec.low), digits = 2)
# 
# icv.high.imput <- signif(subset(out.lst$k4me3$dat.meta, cell == high.var.cell)$cell.var.within.sum.norm, digits = 2)
# icv.low.imput <- signif(subset(out.lst$k4me3$dat.meta, cell == low.var.cell)$cell.var.within.sum.norm, digits = 2)
# 
# ncuts.high <- signif(subset(out.lst$k4me3$dat.meta, cell == high.var.cell)$ncuts, digits = 2)
# ncuts.low <- signif(subset(out.lst$k4me3$dat.meta, cell == low.var.cell)$ncuts, digits = 2)
# 
# cvec <- sapply(rownames(mat.filt), GetChromo)
# # jhash <- hash(unique(cvec), seq(unique(cvec)) %% 2)
# cvec.col <- sapply(cvec, function(x) ifelse(x %in% c("chr11", "chr19"), "blue", "red"))
# 
# 
# par(mfrow=c(2,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
#   plot(xvec.high, main = paste("High Var Cell:", icv.high.imput, icv.high, ncuts.high), col = cvec.col, type = "l")
#   plot(xvec.low, main = paste("Low Var Cell:", icv.low.imput, icv.low, ncuts.low), col = cvec.col, type = "l")
# par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# 
# par(mfrow=c(2,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
#   acf(xvec.high)
#   acf(xvec.low)
# par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# 
# # calculate icv for every cell
# cell.means <- data.frame(cell = cells.keep, jmean = colMeans(mat[, cells.keep]))
# dat.var.raw <- CalculateVarAll(as.matrix(mat[, cells.keep]), "") %>%
#   dplyr::rename(icv.raw = cell.var.within.sum.norm) %>%
# 
# dat.var.merge <- left_join(out.lst$k4me3$dat.meta, subset(dat.var.raw, select = c(cell, icv.raw)))
# dat.var.merge <- left_join(dat.var.merge, cell.means)
# 
# # confounded by cell size?
# ggplot(dat.var.merge, aes(x = icv.raw, y = cell.var.within.sum.norm, color = ncuts)) + geom_point() + scale_color_viridis_c()
# 
# ggplot(dat.var.merge, aes(x = sqrt(icv.raw) / jmean, y = cell.var.within.sum.norm, color = log10(ncuts))) + geom_point() + scale_color_viridis_c()
# 
# 
# 
# # use multinom model? 
# 
# 
# # downsamp? 
# jprop <- min(colSums(mat.filt)) / colSums(mat.filt)
# mat.filt.downsamp <- downsampleMatrix(mat.filt, jprop)
# 
# par(mfrow=c(2,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
#   plot(mat.filt.downsamp[, high.var.cell], main = "High Var Cell")
#   plot(mat.filt.downsamp[, low.var.cell], main = "Low Var Cell")
# par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# 
# # if no reads on chr1, what is variance?
# (jcells.filt <- which(colSums(mat.filt) == 0))


# Check with LSI ----------------------------------------------------------

# 
# 
# mat.filt.allterms <- mat[, cells.keep.high]
# library(irlba)
# mat.filt.allterms.downsamp <- downsampleMatrix(mat.filt.allterms, prop = min(colSums(mat.filt.allterms)) / colSums(mat.filt.allterms))
# # lsi.out <- RunLSI(as.matrix(mat.filt.allterms))
# lsi.out <- RunLSI(as.matrix(mat.filt.allterms.downsamp))
# 
# topics.mat <- lsi.out$u
# jsettings <- umap.defaults
# jsettings$n_neighbors <- 30
# jsettings$min_dist <- 0.1
# jsettings$random_state <- 123
# umap.out <- umap(topics.mat, config = jsettings)
# dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
# dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long) %>%
#   left_join(., out.lst$k4me3$dat.meta)
# cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)
# 
# ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_viridis_c()
# 
# 



# 
# # check downsamp
# 
# inf.ds <- "/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_BinFiltCellFilt2/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.TRUE.no_filt.downsamp_2000/lda_out_meanfilt.B6_H3K4me3_pcutoff_0.CountThres0.K-25_30_35_50.Robj"
# load(inf.ds, v=T)
# 
# out.lda[[1]]@terms <- paste("chr", out.lda[[1]]@terms, sep = "")
# tm.result <- posterior(out.lda[[1]])
# dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
# 
# dat.var <- CalculateVarAll(dat.impute.log = dat.impute.log, jchromos = jchromos)
# 
# topics.mat <- posterior(out.lda[[1]])$topics
# jsettings <- umap.defaults
# jsettings$n_neighbors <- 30
# jsettings$min_dist <- 0.1
# jsettings$random_state <- 123
# umap.out <- umap(topics.mat, config = jsettings)
# dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
# dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long) %>%
#   left_join(., dat.var)
# cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)
# 
# ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
#   geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_fill_viridis_c(direction = -1) 
# 



