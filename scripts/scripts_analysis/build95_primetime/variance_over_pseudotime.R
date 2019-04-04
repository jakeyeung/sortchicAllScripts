# Jake Yeung
# Date of Creation: 2019-04-03
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/variance_over_pseudotime.R
# Pseudotime  =over variance


rm(list=ls())

library(dplyr)
library(ggplot2)
library(tidyr)

source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/FourierFunctions.R")

GetDiffRelToCell <- function(imputed.dat, gstr, trajs, trajname){
  hsc.cell <- (trajs[[jmark]][[trajname]] %>% arrange(lambda) %>% dplyr::top_n(-1))$cell[[1]]
  jsub.all <- MatToLong(imputed.dat, gstr = "", cells.vec=NULL)
  jsub.hsc <- jsub.all %>% filter(cell == hsc.cell)
  # plot by reference to stem cell 
  jsub.hsc.ref <- jsub.hsc %>% rename(exprs.ref = exprs) %>% select(-cell, -start, -end, -pos, -chromo)
  jsub.ref <- left_join(jsub.all, jsub.hsc.ref)
  # do the difference over pseudotime?? 
  jsub.ref.sum <- jsub.ref %>%
    rowwise() %>%
    mutate(exprs.diff = log2(exprs) - log2(exprs.ref)) %>%
    group_by(cell, chromo) %>%
    summarise(exprs.diff.med = median(abs(exprs.diff)))
  # join to UMAP 
  jsub.ref.merge <- left_join(jsub.ref.sum %>% dplyr::select(cell, exprs.diff.med), dat.umap.long)
  return(jsub.ref.merge)
}

GetCellSd <- function(dat.mat, grep.str, log2.scale = TRUE, fn = sd){
  # calculate standard deviation from matrix
  row.filt.indx <- grepl(grep.str, rownames(dat.mat))
  if (log2.scale){
    cell.sd.df <- data.frame(cell = colnames(dat.mat[row.filt.indx, ]),
                             cell.sd = apply(log2(dat.mat[row.filt.indx, ]), 2, fn))
  } else {
    cell.sd.df <- data.frame(cell = colnames(dat.mat[row.filt.indx, ]),
                             cell.sd = apply(dat.mat[row.filt.indx, ], 2, fn))
  }
  cell.sd.df$label <- grep.str
  return(cell.sd.df)
}

MatToLong <- function(imputed.dat, gstr, cells.vec = NULL){
  if (!is.null(cells.vec)){
    jsub <- as.data.frame(imputed.dat[grepl(gstr, rownames(imputed.dat)), cells.vec])
  } else {
    jsub <- as.data.frame(imputed.dat[grepl(gstr, rownames(imputed.dat)), ])
  }
  jsub$coord <- rownames(jsub)
  jsub$start <- as.numeric(sapply(jsub$coord, GetStart))
  jsub$end <- as.numeric(sapply(jsub$coord, GetStart))
  jsub$pos <- jsub$start + (jsub$end - jsub$start) / 2
  jsub$chromo <- sapply(jsub$coord, GetChromo)
  jsub <- gather(jsub, key = "cell", value = "exprs", c(-coord, -start, -end, -pos))
  jsub <- jsub %>% arrange(desc(pos))
  return(jsub)
}

MergeSdWithPseudotime <- function(dat.umap.long.trajs, tm.result.lst, trajs, jmark, jtraj, grep.strs, jscale=TRUE, jfn = mad){
  imputed.dat <- t(tm.result.lst[[jmark]]$terms) %*% t(tm.result.lst[[jmark]]$topics)
  dat.umap.long <- dat.umap.long.trajs[[jmark]]
  cell.sd.df.long <- lapply(grep.strs, function(grep.str){
    return(GetCellSd(imputed.dat, grep.str, log2.scale = jscale, fn = jfn))
  }) %>%
    bind_rows()
  dat.umap.filt <- left_join(dat.umap.long, cell.sd.df.long)
  # add a trajectory
  dat.umap.filt <- left_join(trajs[[jmark]][[jtraj]] %>% dplyr::select(cell, lambda), dat.umap.filt, by = "cell")
  return(dat.umap.filt)
}

# Load data ---------------------------------------------------------------

inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.RData"

load(inf, v=T)


# Load constants ----------------------------------------------------------


jcol <- c("gray80", "gray50", "darkblue")

# Do H3K27me3 -------------------------------------------------------------


jmark <- "H3K27me3"
imputed.dat <- t(tm.result.lst[[jmark]]$terms) %*% t(tm.result.lst[[jmark]]$topics)
dat.umap.long <- dat.umap.long.trajs[[jmark]]

grep.strs <- paste("chr", c(seq(21)), ":", sep = "")

cell.sd.df.long <- lapply(grep.strs, function(grep.str){
  return(GetCellSd(imputed.dat, grep.str, log2.scale = TRUE))
}) %>%
  bind_rows()

dat.umap.filt <- left_join(dat.umap.long, cell.sd.df.long)

m.chr <- ggplot(dat.umap.filt, aes(x = umap1, y = umap2, color = cell.sd)) + 
  geom_point(size = 0.2)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = 1, limits = c(0, 3))

print(m.chr)


# what about genome wide
cell.sd.genomewide <- GetCellSd(imputed.dat, "", log2.scale=TRUE)

dat.umap.filt.gw <- left_join(dat.umap.long, cell.sd.genomewide)

m.gw <- ggplot(dat.umap.filt.gw, aes(x = umap1, y = umap2, color = cell.sd)) + 
  geom_point(size = 2)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = 1.5) + ggtitle("Genome wide")

print(m.gw)

hsc.cell <- (trajs[[jmark]]$myeloid %>% arrange(lambda) %>% dplyr::top_n(-1))$cell[[1]]
myeloid.cell <- (trajs[[jmark]]$myeloid %>% arrange(lambda) %>% dplyr::top_n(1))$cell[[1]]


gstr <- paste0("chr15:")
jsub.all <- as.data.frame(imputed.dat[grepl(gstr, rownames(imputed.dat)), ])


jsub <- as.data.frame(imputed.dat[grepl(gstr, rownames(imputed.dat)), c(hsc.cell, myeloid.cell)])
jsub$coord <- rownames(jsub)
jsub <- gather(jsub, key = "cell", value = "exprs", -coord)
jsub$start <- as.numeric(sapply(jsub$coord, GetStart))
jsub$end <- as.numeric(sapply(jsub$coord, GetStart))
jsub$pos <- jsub$start + (jsub$end - jsub$start) / 2
jsub <- jsub %>% arrange(desc(pos))


ggplot(jsub, aes(x = pos / 10^6, y = log2(exprs))) + 
  geom_point() + 
  facet_wrap(~cell) + 
  ggtitle(gstr) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jsub, aes(x = pos / 10^6, y = log2(exprs), group = cell, color = cell)) + 
  geom_point() + 
  ggtitle(gstr) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jsub %>% group_by(pos) %>% summarise(exprs = diff(log2(exprs))), aes(x = pos / 10^6, y = exprs)) + 
  geom_point() + 
  ggtitle(gstr) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# spatial pattern?

jsub.hsc <- jsub %>% filter(cell == hsc.cell)
jsub.myeloid <- jsub %>% filter(cell == myeloid.cell)
acf(log2(jsub.hsc$exprs), type = "partial")
acf(log2(jsub.hsc$exprs), type = "partial", lag.max = nrow(jsub.hsc))
acf(log2(jsub.myeloid$exprs), type = "partial")
acf(log2(jsub.myeloid$exprs), type = "partial", lag.max = nrow(jsub.myeloid))


m1 <- ggplot(jsub.hsc, aes(x = pos, y = log2(exprs))) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
m2 <- ggplot(jsub.myeloid, aes(x = pos, y = log2(exprs))) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
multiplot(m1 + ggtitle(paste(gstr, jmark)), m2, cols = 2)



head(jsub)


jsub.ref.merge <- GetDiffRelToCell(imputed.dat, gstr = "chr15:", trajs, trajname = "myeloid")

ggplot(jsub.ref.merge, aes(x = umap1, y = umap2, color = exprs.diff.med)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = 0.5) + facet_wrap(~chromo)


# Do H3K9me3  -------------------------------------------------------------


jmark <- "H3K9me3"
imputed.dat <- t(tm.result.lst[[jmark]]$terms) %*% t(tm.result.lst[[jmark]]$topics)
dat.umap.long <- dat.umap.long.trajs[[jmark]]

grep.strs <- paste("chr", c(seq(21)), ":", sep = "")

cell.sd.df.long <- lapply(grep.strs, function(grep.str){
  return(GetCellSd(imputed.dat, grep.str, log2.scale = TRUE))
}) %>%
  bind_rows()
dat.umap.filt <- left_join(dat.umap.long, cell.sd.df.long)

m.chr <- ggplot(dat.umap.filt, aes(x = umap1, y = umap2, color = cell.sd)) + 
  geom_point(size = 0.2)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = 1)
print(m.chr)


# what about genome wide
cell.sd.genomewide <- GetCellSd(imputed.dat, "", log2.scale=TRUE)

dat.umap.filt.gw <- left_join(dat.umap.long, cell.sd.genomewide)

m.gw <- ggplot(dat.umap.filt.gw, aes(x = umap1, y = umap2, color = cell.sd)) + 
  geom_point(size = 2)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = 1.5) + ggtitle("Genome wide")

print(m.gw)



# ggplot(dat.umap.filt %>% filter(label == "chr1:"), aes(x = ))

# do myeloid trajectory

hsc.cell <- (trajs[[jmark]]$myeloid %>% arrange(lambda) %>% dplyr::top_n(-1))$cell[[1]]
myeloid.cell <- (trajs[[jmark]]$myeloid %>% arrange(lambda) %>% dplyr::top_n(1))$cell[[1]]

gstr <- paste0("chr15:")
jsub <- as.data.frame(imputed.dat[grepl(gstr, rownames(imputed.dat)), c(hsc.cell, myeloid.cell)])
jsub$coord <- rownames(jsub)
jsub <- gather(jsub, key = "cell", value = "exprs", -coord)
jsub$start <- as.numeric(sapply(jsub$coord, GetStart))
jsub$end <- as.numeric(sapply(jsub$coord, GetStart))
jsub$pos <- jsub$start + (jsub$end - jsub$start) / 2
jsub <- jsub %>% arrange(desc(pos))

ggplot(jsub, aes(x = pos / 10^6, y = log2(exprs))) + 
  geom_point() + 
  facet_wrap(~cell) + 
  ggtitle(gstr) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


jsub.hsc <- jsub %>% filter(cell == hsc.cell)
jsub.myeloid <- jsub %>% filter(cell == myeloid.cell)
acf(log2(jsub.hsc$exprs), type = "partial")
acf(log2(jsub.myeloid$exprs), type = "partial")

m1 <- ggplot(jsub.hsc, aes(x = pos, y = exprs)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
m2 <- ggplot(jsub.myeloid, aes(x = pos, y = exprs)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
multiplot(m1 + ggtitle(paste(gstr, jmark)), m2, cols = 2)


head(jsub)

subset(jsub, coord == "chr15:103940000-104040000")

jsub.diff <- jsub %>%
  group_by(coord) %>%
  summarise(exprs.diff = abs(diff(log2(exprs)))) %>%
  ungroup() %>%
  summarise(exprs.diff.sum = median(exprs.diff))

print(jsub.diff)


# Do H3K4me1 --------------------------------------------------------------


jmark <- "H3K4me1"
imputed.dat <- t(tm.result.lst[[jmark]]$terms) %*% t(tm.result.lst[[jmark]]$topics)
dat.umap.long <- dat.umap.long.trajs[[jmark]]

grep.strs <- paste("chr", c(seq(21)), ":", sep = "")

jscale <- TRUE
cell.sd.df.long <- lapply(grep.strs, function(grep.str){
  return(GetCellSd(imputed.dat, grep.str, log2.scale = jscale, fn = mad))
}) %>%
  bind_rows()
dat.umap.filt <- left_join(dat.umap.long, cell.sd.df.long)

if (jscale){
  m.chr <- ggplot(dat.umap.filt, aes(x = umap1, y = umap2, color = cell.sd)) + 
    geom_point(size = 0.2)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = 1.5, limits = c(0, 3))
} else {
  m.chr <- ggplot(dat.umap.filt, aes(x = umap1, y = umap2, color = cell.sd)) + 
    geom_point(size = 0.2)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = 5e-6)
}
print(m.chr)


# what about genome wide
cell.sd.genomewide <- GetCellSd(imputed.dat, "", log2.scale=jscale, fn = mad)

dat.umap.filt.gw <- left_join(dat.umap.long, cell.sd.genomewide)

if (jscale){
  m.gw <- ggplot(dat.umap.filt.gw, aes(x = umap1, y = umap2, color = cell.sd)) + 
    geom_point(size = 2)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = 1.5, limits = c(0, 3)) + 
    ggtitle("Genome wide")
} else {
  m.gw <- ggplot(dat.umap.filt.gw, aes(x = umap1, y = umap2, color = cell.sd)) + 
    geom_point(size = 2)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = 5e-6) + 
    ggtitle("Genome wide")
}
print(m.gw)


hsc.cell <- (trajs[[jmark]]$myeloid %>% arrange(lambda) %>% dplyr::top_n(-1))$cell[[1]]
myeloid.cell <- (trajs[[jmark]]$myeloid %>% arrange(lambda) %>% dplyr::top_n(1))$cell[[1]]

gstr <- paste0("chr15:")
jsub <- as.data.frame(imputed.dat[grepl(gstr, rownames(imputed.dat)), c(hsc.cell, myeloid.cell)])
jsub$coord <- rownames(jsub)
jsub <- gather(jsub, key = "cell", value = "exprs", -coord)
jsub$start <- as.numeric(sapply(jsub$coord, GetStart))
jsub$end <- as.numeric(sapply(jsub$coord, GetStart))
jsub$pos <- jsub$start + (jsub$end - jsub$start) / 2
jsub <- jsub %>% arrange(desc(pos))

ggplot(jsub, aes(x = pos / 10^6, 
                 y = exprs)) + 
  geom_point() + 
  facet_wrap(~cell) + 
  ggtitle(gstr) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jsub, aes(x = pos / 10^6, y = log2(exprs), group = cell, color = cell)) + 
  geom_point() + 
  ggtitle(gstr) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jsub %>% group_by(pos) %>% summarise(exprs = diff(log2(exprs))), aes(x = pos / 10^6, y = exprs)) + 
  geom_point() + 
  ggtitle(gstr) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jsub.hsc <- jsub %>% filter(cell == hsc.cell)
jsub.myeloid <- jsub %>% filter(cell == myeloid.cell)
acf(log2(jsub.hsc$exprs), type = "partial")
acf(log2(jsub.myeloid$exprs), type = "partial")

m1 <- ggplot(jsub.hsc, aes(x = pos, y = exprs)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
m2 <- ggplot(jsub.myeloid, aes(x = pos, y = exprs)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
multiplot(m1 + ggtitle(paste(gstr, jmark)), m2, cols = 2)

head(jsub)

subset(jsub, coord == "chr15:103940000-104040000")

jsub.diff <- jsub %>%
  group_by(coord) %>%
  summarise(exprs.diff = abs(diff(log2(exprs)))) %>%
  ungroup() %>%
  summarise(exprs.diff.sum = median(exprs.diff))

print(jsub.diff)

# Do H3K4me3 --------------------------------------------------------------

jmark <- "H3K4me3"
imputed.dat <- t(tm.result.lst[[jmark]]$terms) %*% t(tm.result.lst[[jmark]]$topics)
dat.umap.long <- dat.umap.long.trajs[[jmark]]

grep.strs <- paste("chr", c(seq(21)), ":", sep = "")

cell.sd.df.long <- lapply(grep.strs, function(grep.str){
  return(GetCellSd(imputed.dat, grep.str, log2.scale = TRUE))
}) %>%
  bind_rows()
dat.umap.filt <- left_join(dat.umap.long, cell.sd.df.long)

m.chr <- ggplot(dat.umap.filt, aes(x = umap1, y = umap2, color = cell.sd)) + 
  geom_point(size = 0.2)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = 1)
print(m.chr)

# what about genome wide
cell.sd.genomewide <- GetCellSd(imputed.dat, "", log2.scale=TRUE)

dat.umap.filt.gw <- left_join(dat.umap.long, cell.sd.genomewide)

m.gw <- ggplot(dat.umap.filt.gw, aes(x = umap1, y = umap2, color = cell.sd)) + 
  geom_point(size = 2)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = 1.5) + ggtitle("Genome wide")

print(m.gw)



hsc.cell <- (trajs[[jmark]]$myeloid %>% arrange(lambda) %>% dplyr::top_n(-1))$cell[[1]]
myeloid.cell <- (trajs[[jmark]]$myeloid %>% arrange(lambda) %>% dplyr::top_n(1))$cell[[1]]

gstr <- paste0("chr15:")
jsub <- as.data.frame(imputed.dat[grepl(gstr, rownames(imputed.dat)), c(hsc.cell, myeloid.cell)])
jsub$coord <- rownames(jsub)
jsub <- gather(jsub, key = "cell", value = "exprs", -coord)
jsub$start <- as.numeric(sapply(jsub$coord, GetStart))
jsub$end <- as.numeric(sapply(jsub$coord, GetStart))
jsub$pos <- jsub$start + (jsub$end - jsub$start) / 2
jsub <- jsub %>% arrange(desc(pos))

ggplot(jsub, aes(x = pos / 10^6, y = log2(exprs))) + 
  geom_point() + 
  facet_wrap(~cell) + 
  ggtitle(gstr) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(jsub, aes(x = pos / 10^6, y = log2(exprs), group = cell, color = cell)) + 
  geom_point() + 
  ggtitle(gstr) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jsub %>% group_by(pos) %>% summarise(exprs = diff(log2(exprs))), aes(x = pos / 10^6, y = exprs)) + 
  geom_point() + 
  ggtitle(gstr) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# spatial pattern?

jsub.hsc <- jsub %>% filter(cell == hsc.cell)
jsub.myeloid <- jsub %>% filter(cell == myeloid.cell)
jout <- acf(log2(jsub.hsc$exprs), type = "correlation", lag.max = nrow(jsub.hsc))
# jout <- acf(log2(jsub.hsc$exprs), type = "correlation", lag.max = 5, plot = TRUE)
acf(log2(jsub.hsc$exprs), type = "partial", lag.max = nrow(jsub.hsc))
acf(log2(jsub.myeloid$exprs), type = "correlation", lag.max = nrow(jsub.myeloid))
acf(log2(jsub.myeloid$exprs), type = "partial", lag.max = nrow(jsub.myeloid))


m1 <- ggplot(jsub.hsc, aes(x = pos, y = log2(exprs))) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
m2 <- ggplot(jsub.myeloid, aes(x = pos, y = log2(exprs))) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# m1 <- ggplot(jsub.hsc, aes(x = rank(pos), y = exprs)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# m2 <- ggplot(jsub.myeloid, aes(x = rank(pos), y = exprs)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
multiplot(m1 + ggtitle(paste(gstr, jmark)), m2, cols = 2)



# plot exprs vs exprslag1

# jlag <- 2000
jlag <- 24
jlag <- 1010 - 1
jsub.hsc$exprslag1 <- c(rep(NA, jlag), jsub.hsc$exprs[1:(length(jsub.hsc$exprs) - jlag)])
ggplot(jsub.hsc, aes(x = exprs, exprslag1)) + geom_point(alpha = 0.2) + theme_bw()

# ggplot(jsub.myeloid, aes(x = exprs, exprslag1)) + geom_point(alpha = 0.2) + theme_bw()

# ar.fit <- ar(x = scale(jsub$exprs), order.max = 3)
afits.myeloid <- lapply(seq(2, 12), function(x) arima(x = jsub.myeloid$exprs, order = c(x,0,0)))
afits.hsc <- lapply(seq(2, 12), function(x) arima(x = jsub.hsc$exprs, order = c(x,0,0)))

bics.hsc <- lapply(afits.hsc, BIC)
bics.myeloid <- lapply(afits.myeloid, BIC)

# spatial pattern by periodogram
samp.interval <- 20000  # 20kb
pout <- CalculatePeriodogram(scale(jsub.hsc$exprs))
max.freqs <- FindMaxFreqs(pout$freq, pout$p.scaled)
Tmax <- samp.interval / max.freqs
PlotPeriodogram(Frequency = pout$freq, Periodogram = pout$p.scaled)

plot(x = (samp.interval / pout$freq),  y = pout$p.scaled, type = "o", log = "x")


# compare HSC state vs differentiated state
head(jsub)

subset(jsub, coord == "chr15:103940000-104040000")

jsub.diff <- jsub %>%
  group_by(coord) %>%
  summarise(exprs.diff = abs(diff(log2(exprs)))) %>%
  ungroup() %>%
  summarise(exprs.diff.sum = median(exprs.diff))


# Compare MADs in repressive and active marks -----------------------------





# jchromo <- "chr15:"
jtraj <- "myeloid"
jmark <- "H3K4me3"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")

grep.strs <- paste("chr", c(seq(21)), ":", sep = "")


grep.strs <- ""

dat.umap.filt.lst <- lapply(jmarks, function(jmark) MergeSdWithPseudotime(dat.umap.long.trajs, tm.result.lst, trajs, 
                                                                          jmark, jtraj, grep.strs, jscale=TRUE, jfn=mad))

dat.umap.filt.long <- dat.umap.filt.lst %>% bind_rows()
dat.umap.filt.long$mark <- factor(as.character(dat.umap.filt.long$mark), levels = jmarks)
# plot everything together

ggplot(dat.umap.filt.long, aes(x = lambda, y = cell.sd, group = mark, color = mark)) + 
  geom_point(size = 1, alpha = 1) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~label)

# plot one chromosome across marks 

# jchromo <- "chr19:"
jchromo <- ""
ggplot(dat.umap.filt.long %>% filter(label == jchromo), aes(x = lambda, y = cell.sd, group = mark, color = mark)) + 
  geom_point(size = 0.5) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark) + ggtitle(jchromo)


# Play with difference from stem cell state -------------------------------



# Test noise --------------------------------------------------------------


# square wave can induce temporal correlations ?
library(zoo)
x <- rnorm(n = 5000, mean = 0, sd = 2)
x.smooth <- rollapply(x, width = 5, by = 1, FUN = mean, align = "left")

plot(x)
acf(x, type = "correlation")
acf(x, type = "partial")
plot(x.smooth)
acf(x.smooth, type = "correlation")
acf(x.smooth, type = "partial")

