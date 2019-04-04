# Jake Yeung
# Date of Creation: 2019-04-03
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/variance_over_pseudotime_plots.R
# Plot the outputs 

rm(list=ls())

tstart <- Sys.time()

library(dplyr)
library(ggplot2)
library(tidyr)
library(JFuncs)

source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/FourierFunctions.R")

GetDiffRelToCell <- function(imputed.dat, gstr, trajs, trajname, dat.umap.long){
  hsc.cell <- (trajs[[jmark]][[trajname]] %>% arrange(lambda) %>% dplyr::top_n(-1))$cell[[1]]
  jsub.all <- MatToLong(imputed.dat, gstr = gstr, cells.vec=NULL)
  jsub.hsc <- jsub.all %>% filter(cell == hsc.cell)
  # plot by reference to stem cell 
  jsub.hsc.ref <- jsub.hsc %>% rename(exprs.ref = exprs) %>% select(-cell, -start, -end, -pos, -chromo)
  jsub.ref <- left_join(jsub.all, jsub.hsc.ref)
  # do the difference over pseudotime?? 
  jsub.ref$exprs.diff <- log2(jsub.ref$exprs) - log2(jsub.ref$exprs.ref)
  jsub.ref.sum <- jsub.ref %>%
    group_by(cell) %>%
    summarise(exprs.diff.med = median(abs(exprs.diff)))
  # join to UMAP 
  jsub.ref.merge <- left_join(jsub.ref.sum %>% dplyr::select(cell, exprs.diff.med), dat.umap.long) %>%
    mutate(label = gstr)
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
  jsub <- gather(jsub, key = "cell", value = "exprs", c(-coord, -start, -end, -pos, -chromo))
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


CalculateACF <- function(jsub.hsc, jstep = 20000, jtype = "correlation", jmain = "Title", show.plot = TRUE){
  # impute missing positions with minimum value
  # impute missing bins with minimum value
  jcells <- unique(jsub.hsc$cell)
  pos.start <- min(jsub.hsc$pos)
  pos.end <- max(jsub.hsc$pos)
  # jstep <- 20000
  pos.vec <- seq(pos.start, pos.end, jstep)
  jsub.impute.vec <- data.frame(pos = rep(pos.vec, length(jcells)), cell = rep(jcells, each = length(pos.vec)))
  jsub.impute.vec <- left_join(jsub.impute.vec, jsub.hsc %>% dplyr::select(c(chromo, pos, cell, exprs)))
  # jsub.impute.vec$exprs[which(is.na(jsub.impute.vec$exprs))] <- min(jsub.hsc$exprs)
  
  acf.out <- acf(log2(jsub.impute.vec$exprs), type = jtype, lag.max = nrow(jsub.impute.vec), na.action = na.pass, main = jmain, plot = show.plot)
  acf.out$lag.stepadj <- acf.out$lag * jstep
  return(acf.out)
}


# Load constants ----------------------------------------------------------


jcol <- c("gray80", "gray50", "darkblue")
grep.strs <- paste("chr", c(seq(21)), ":", sep = "")

jalpha <- 0.5

pseudo <- 0
jscale <- 1

mdpt.sd <- 1
lims.sd <- c(0, 3)
mdpt.fc <- 0.75
lims.fc <- c(0, 3)

jsize.facet <- 0.2
gw.jsize.facet <- 2

do.plots <- FALSE

jstep <- 20000
jtype <- "correlation"


# Load data ---------------------------------------------------------------

inf <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.RData"

load(inf, v=T)


# Do H3K27me3 -------------------------------------------------------------

jmark <- "H3K27me3"
trajname <- "myeloid"

print(paste(jmark, trajname))

imputed.dat <- t(tm.result.lst[[jmark]]$terms) %*% t(tm.result.lst[[jmark]]$topics)
dat.umap.long <- dat.umap.long.trajs[[jmark]]


if (do.plots){
  pdf("~/data/scchic/pdfs/variance_over_pseudotime_plots.pdf", useDingbats = FALSE)
}

cell.sd.df.long <- lapply(grep.strs, function(grep.str){
  return(GetCellSd(jscale * (imputed.dat + pseudo), grep.str, log2.scale = TRUE, fn = sd))
}) %>%
  bind_rows()

dat.umap.filt <- left_join(dat.umap.long, cell.sd.df.long)

print(range(dat.umap.filt$cell.sd))

m.chr <- ggplot(dat.umap.filt, aes(x = umap1, y = umap2, color = cell.sd)) + 
  geom_point(size = jsize.facet)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.sd, limits = lims.sd) + 
  ggtitle(jmark, paste0(deparse(substitute(sd)), " across chromosome"))
print(m.chr)

# what about genome wide
cell.sd.genomewide <- GetCellSd(imputed.dat, "", log2.scale=TRUE)
dat.umap.filt.gw <- left_join(dat.umap.long, cell.sd.genomewide)

m.gw <- ggplot(dat.umap.filt.gw, aes(x = umap1, y = umap2, color = cell.sd)) + 
  geom_point(size = gw.jsize.facet)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.sd, limits = lims.sd) + 
  ggtitle(jmark, paste0(deparse(substitute(sd)), " genome wide"))

print(m.gw)

  
  # Highlight differences for two representative cells on a representative chromosome 
  hsc.cell <- (trajs[[jmark]][[trajname]] %>% arrange(lambda) %>% dplyr::top_n(-1))$cell[[1]]
  diff.cell <- (trajs[[jmark]][[trajname]] %>% arrange(lambda) %>% dplyr::top_n(1))$cell[[1]]
  
  gstr <- paste0("chr15:")
  jsub <- MatToLong(imputed.dat, gstr, cells.vec = c(hsc.cell, diff.cell))
  m.spatial <- ggplot(jsub, aes(x = pos / 10^6, y = log2(exprs))) + 
    geom_line(alpha = jalpha) + 
    facet_wrap(~cell) + 
    ggtitle(paste(jmark, gstr)) + 
    xlab("MB") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  m.spatial.merged <- ggplot(jsub, aes(x = pos / 10^6, y = log2(exprs), group = cell, color = cell)) + 
    geom_line(alpha = jalpha) + 
    ggtitle(paste(jmark, gstr)) + 
    xlab("MB") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  m.spatial.log2fc <- ggplot(jsub %>% group_by(pos) %>% summarise(exprs = diff(log2(exprs))), aes(x = pos / 10^6, y = exprs)) + 
    geom_line(alpha = jalpha) + 
    ggtitle(paste(jmark, gstr)) + 
    xlab("MB") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ylab("log2 Fold Change")
  
  print(m.spatial)
  print(m.spatial.merged)
  print(m.spatial.log2fc)
  
  
  # spatial pattern?
  
  jsub.hsc <- jsub %>% filter(cell == hsc.cell)
  jsub.myeloid <- jsub %>% filter(cell == diff.cell)
  acf(log2(jsub.hsc$exprs), type = "partial", main = paste(jmark, gstr, "HSC Cell autocorrelation"))
  acf(log2(jsub.hsc$exprs), type = "partial", lag.max = nrow(jsub.hsc), main = paste(jmark, gstr, "HSC Cell autocorrelation"))
  acf(log2(jsub.myeloid$exprs), type = "partial", main = paste(jmark, gstr, "Myeloid cell autocorrelation"))
  acf(log2(jsub.myeloid$exprs), type = "partial", lag.max = nrow(jsub.myeloid), main = paste(jmark, gstr, "Myeloid cell autocorrelation"))
  
  
  # Plot the median log2 fold change relative to HSC cell: for one chromo
  
  jsub.ref.merge <- lapply(grep.strs, function(gstr) GetDiffRelToCell(imputed.dat, gstr = gstr, trajs, trajname = trajname, dat.umap.long)) %>%
    bind_rows() 
  m.mad <- ggplot(jsub.ref.merge, aes(x = umap1, y = umap2, color = exprs.diff.med)) + geom_point(size = jsize.facet) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.fc, limits = lims.fc) + facet_wrap(~label) + ggtitle(jmark, "Median log2 fold change with HSC")
  print(m.mad)
  
  # do genome-wide?
  jsub.ref.merge.gw <- lapply(c(""), function(gstr) GetDiffRelToCell(imputed.dat, gstr = gstr, trajs, trajname = trajname, dat.umap.long)) %>%
    bind_rows() 
  m.mad.gw <- ggplot(jsub.ref.merge.gw, aes(x = umap1, y = umap2, color = exprs.diff.med)) + geom_point(size = gw.jsize.facet) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.fc, limits = lims.fc) + facet_wrap(~label) + ggtitle(jmark, "Median log2 fold change with HSC")
  print(m.mad.gw)
  
  print(range(jsub.ref.merge.gw$exprs.diff.med))
  # plot along pseudotime? 
  traj.sub <- trajs[[jmark]][[trajname]]
  # add exprs.diff.med
  traj.sub <- left_join(traj.sub, jsub.ref.merge.gw %>% dplyr::select(cell, exprs.diff.med), by = c("cell"))
  
  m.mad.traj <- ggplot(traj.sub, aes(x = lambda, y = exprs.diff.med)) + geom_point(alpha = 0.1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlab("Pseudotime") + ylab("Median Log2 FC from Prog Cell") + 
    ggtitle(jmark, paste(trajname, "Genome-wide"))
  print(m.mad.traj)
  

# H3K9me3 -----------------------------------------------------------------


jmark <- "H3K9me3"
trajname <- "myeloid"

print(paste(jmark, trajname))

imputed.dat <- t(tm.result.lst[[jmark]]$terms) %*% t(tm.result.lst[[jmark]]$topics)
dat.umap.long <- dat.umap.long.trajs[[jmark]]



cell.sd.df.long <- lapply(grep.strs, function(grep.str){
  return(GetCellSd(jscale * (imputed.dat + pseudo), grep.str, log2.scale = TRUE, fn = sd))
}) %>%
  bind_rows()

dat.umap.filt <- left_join(dat.umap.long, cell.sd.df.long)

print(range(dat.umap.filt$cell.sd))

m.chr <- ggplot(dat.umap.filt, aes(x = umap1, y = umap2, color = cell.sd)) + 
  geom_point(size = jsize.facet)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.sd, limits = lims.sd) + 
  ggtitle(jmark, paste0(deparse(substitute(sd)), " across chromosome"))
print(m.chr)

# what about genome wide
cell.sd.genomewide <- GetCellSd(imputed.dat, "", log2.scale=TRUE)
dat.umap.filt.gw <- left_join(dat.umap.long, cell.sd.genomewide)

m.gw <- ggplot(dat.umap.filt.gw, aes(x = umap1, y = umap2, color = cell.sd)) + 
  geom_point(size = gw.jsize.facet)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.sd, limits = lims.sd) + 
  ggtitle(jmark, paste0(deparse(substitute(sd)), " genome wide"))

print(m.gw)


# Highlight differences for two representative cells on a representative chromosome 
hsc.cell <- (trajs[[jmark]][[trajname]] %>% arrange(lambda) %>% dplyr::top_n(-1))$cell[[1]]
diff.cell <- (trajs[[jmark]][[trajname]] %>% arrange(lambda) %>% dplyr::top_n(1))$cell[[1]]

gstr <- paste0("chr15:")
jsub <- MatToLong(imputed.dat, gstr, cells.vec = c(hsc.cell, diff.cell))
m.spatial <- ggplot(jsub, aes(x = pos / 10^6, y = log2(exprs))) + 
  geom_line(alpha = jalpha) + 
  facet_wrap(~cell) + 
  ggtitle(paste(jmark, gstr)) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
m.spatial.merged <- ggplot(jsub, aes(x = pos / 10^6, y = log2(exprs), group = cell, color = cell)) + 
  geom_line(alpha = jalpha) + 
  ggtitle(paste(jmark, gstr)) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
m.spatial.log2fc <- ggplot(jsub %>% group_by(pos) %>% summarise(exprs = diff(log2(exprs))), aes(x = pos / 10^6, y = exprs)) + 
  geom_line(alpha = jalpha) + 
  ggtitle(paste(jmark, gstr)) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ylab("log2 Fold Change")

print(m.spatial)
print(m.spatial.merged)
print(m.spatial.log2fc)


# spatial pattern?

jsub.hsc <- jsub %>% filter(cell == hsc.cell)
jsub.myeloid <- jsub %>% filter(cell == diff.cell)
acf(log2(jsub.hsc$exprs), type = "partial", main = paste(jmark, gstr, "HSC Cell autocorrelation"))
acf(log2(jsub.hsc$exprs), type = "partial", lag.max = nrow(jsub.hsc), main = paste(jmark, gstr, "HSC Cell autocorrelation"))
acf(log2(jsub.myeloid$exprs), type = "partial", main = paste(jmark, gstr, "Myeloid cell autocorrelation"))
acf(log2(jsub.myeloid$exprs), type = "partial", lag.max = nrow(jsub.myeloid), main = paste(jmark, gstr, "Myeloid cell autocorrelation"))


# Plot the median log2 fold change relative to HSC cell: for one chromo

jsub.ref.merge <- lapply(grep.strs, function(gstr) GetDiffRelToCell(imputed.dat, gstr = gstr, trajs, trajname = trajname, dat.umap.long)) %>%
  bind_rows() 
m.mad <- ggplot(jsub.ref.merge, aes(x = umap1, y = umap2, color = exprs.diff.med)) + geom_point(size = jsize.facet) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.fc, limits = lims.fc) + facet_wrap(~label) + ggtitle(jmark, "Median log2 fold change with HSC")
print(m.mad)

# do genome-wide?
jsub.ref.merge.gw <- lapply(c(""), function(gstr) GetDiffRelToCell(imputed.dat, gstr = gstr, trajs, trajname = trajname, dat.umap.long)) %>%
  bind_rows() 
m.mad.gw <- ggplot(jsub.ref.merge.gw, aes(x = umap1, y = umap2, color = exprs.diff.med)) + geom_point(size = gw.jsize.facet) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.fc, limits = lims.fc) + facet_wrap(~label) + ggtitle(jmark, "Median log2 fold change with HSC")
print(m.mad.gw)

print(range(jsub.ref.merge.gw$exprs.diff.med))
# plot along pseudotime? 
traj.sub <- trajs[[jmark]][[trajname]]
# add exprs.diff.med
traj.sub <- left_join(traj.sub, jsub.ref.merge.gw %>% dplyr::select(cell, exprs.diff.med), by = c("cell"))

m.mad.traj <- ggplot(traj.sub, aes(x = lambda, y = exprs.diff.med)) + geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Pseudotime") + ylab("Median Log2 FC from Prog Cell") + 
  ggtitle(jmark, paste(trajname, "Genome-wide"))
print(m.mad.traj)


# H3K4me1 -----------------------------------------------------------------


jmark <- "H3K4me1"
trajname <- "bcell"

print(paste(jmark, trajname))

imputed.dat <- t(tm.result.lst[[jmark]]$terms) %*% t(tm.result.lst[[jmark]]$topics)
dat.umap.long <- dat.umap.long.trajs[[jmark]]



cell.sd.df.long <- lapply(grep.strs, function(grep.str){
  return(GetCellSd(jscale * (imputed.dat + pseudo), grep.str, log2.scale = TRUE, fn = sd))
}) %>%
  bind_rows()

dat.umap.filt <- left_join(dat.umap.long, cell.sd.df.long)

print(range(dat.umap.filt$cell.sd))

m.chr <- ggplot(dat.umap.filt, aes(x = umap1, y = umap2, color = cell.sd)) + 
  geom_point(size = jsize.facet)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.sd, limits = lims.sd) + 
  ggtitle(jmark, paste0(deparse(substitute(sd)), " across chromosome"))
print(m.chr)

# what about genome wide
cell.sd.genomewide <- GetCellSd(imputed.dat, "", log2.scale=TRUE)
dat.umap.filt.gw <- left_join(dat.umap.long, cell.sd.genomewide)

m.gw <- ggplot(dat.umap.filt.gw, aes(x = umap1, y = umap2, color = cell.sd)) + 
  geom_point(size = gw.jsize.facet)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.sd, limits = lims.sd) + 
  ggtitle(jmark, paste0(deparse(substitute(sd)), " genome wide"))

print(m.gw)


# Highlight differences for two representative cells on a representative chromosome 
hsc.cell <- (trajs[[jmark]][[trajname]] %>% arrange(lambda) %>% dplyr::top_n(-1))$cell[[1]]
diff.cell <- (trajs[[jmark]][[trajname]] %>% arrange(lambda) %>% dplyr::top_n(1))$cell[[1]]

gstr <- paste0("chr15:")
jsub <- MatToLong(imputed.dat, gstr, cells.vec = c(hsc.cell, diff.cell))
m.spatial <- ggplot(jsub, aes(x = pos / 10^6, y = log2(exprs))) + 
  geom_line(alpha = jalpha) + 
  facet_wrap(~cell) + 
  ggtitle(paste(jmark, gstr)) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
m.spatial.merged <- ggplot(jsub, aes(x = pos / 10^6, y = log2(exprs), group = cell, color = cell)) + 
  geom_line(alpha = jalpha) + 
  ggtitle(paste(jmark, gstr)) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
m.spatial.log2fc <- ggplot(jsub %>% group_by(pos) %>% summarise(exprs = diff(log2(exprs))), aes(x = pos / 10^6, y = exprs)) + 
  geom_line(alpha = jalpha) + 
  ggtitle(paste(jmark, gstr)) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ylab("log2 Fold Change")

print(m.spatial)
print(m.spatial.merged)
print(m.spatial.log2fc)


# spatial pattern?

jsub.hsc <- jsub %>% filter(cell == hsc.cell)
jsub.myeloid <- jsub %>% filter(cell == diff.cell)
acf(log2(jsub.hsc$exprs), type = "partial", main = paste(jmark, gstr, "HSC Cell autocorrelation"))
acf(log2(jsub.hsc$exprs), type = "partial", lag.max = nrow(jsub.hsc), main = paste(jmark, gstr, "HSC Cell autocorrelation"))
acf(log2(jsub.myeloid$exprs), type = "partial", main = paste(jmark, gstr, "Myeloid cell autocorrelation"))
acf(log2(jsub.myeloid$exprs), type = "partial", lag.max = nrow(jsub.myeloid), main = paste(jmark, gstr, "Myeloid cell autocorrelation"))


# Plot the median log2 fold change relative to HSC cell: for one chromo

jsub.ref.merge <- lapply(grep.strs, function(gstr) GetDiffRelToCell(imputed.dat, gstr = gstr, trajs, trajname = trajname, dat.umap.long)) %>%
  bind_rows() 
m.mad <- ggplot(jsub.ref.merge, aes(x = umap1, y = umap2, color = exprs.diff.med)) + geom_point(size = jsize.facet) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.fc, limits = lims.fc) + facet_wrap(~label) + ggtitle(jmark, "Median log2 fold change with HSC")
print(m.mad)


# do genome-wide?
jsub.ref.merge.gw <- lapply(c(""), function(gstr) GetDiffRelToCell(imputed.dat, gstr = gstr, trajs, trajname = trajname, dat.umap.long)) %>%
  bind_rows() 
m.mad.gw <- ggplot(jsub.ref.merge.gw, aes(x = umap1, y = umap2, color = exprs.diff.med)) + geom_point(size = gw.jsize.facet) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.fc, limits = lims.fc) + facet_wrap(~label) + ggtitle(jmark, "Median log2 fold change with HSC")
print(m.mad.gw)

print(range(jsub.ref.merge.gw$exprs.diff.med))
# plot along pseudotime? 
traj.sub <- trajs[[jmark]][[trajname]]
# add exprs.diff.med
traj.sub <- left_join(traj.sub, jsub.ref.merge.gw %>% dplyr::select(cell, exprs.diff.med), by = c("cell"))

m.mad.traj <- ggplot(traj.sub, aes(x = lambda, y = exprs.diff.med)) + geom_point(alpha = 0.1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Pseudotime") + ylab("Median Log2 FC from Prog Cell") + 
  ggtitle(jmark, "Genome-wide")
print(m.mad.traj)

ggplot(traj.sub, aes(x = umap1, y = umap2, color = exprs.diff.med)) + geom_point(alpha = 0.1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # xlab("Pseudotime") + ylab("Median Log2 FC from Prog Cell") + 
  ggtitle(jmark, "Genome-wide")


# do genome-wide?
jsub.ref.merge.gw <- lapply(c(""), function(gstr) GetDiffRelToCell(imputed.dat, gstr = gstr, trajs, trajname = trajname, dat.umap.long)) %>%
  bind_rows() 
m.mad.gw <- ggplot(jsub.ref.merge.gw, aes(x = umap1, y = umap2, color = exprs.diff.med)) + geom_point(size = gw.jsize.facet) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.fc, limits = lims.fc) + facet_wrap(~label) + ggtitle(jmark, "Median log2 fold change with HSC")
print(m.mad.gw)

print(range(jsub.ref.merge.gw$exprs.diff.med))
# plot along pseudotime? 
traj.sub <- trajs[[jmark]][[trajname]]
# add exprs.diff.med
traj.sub <- left_join(traj.sub, jsub.ref.merge.gw %>% dplyr::select(cell, exprs.diff.med), by = c("cell"))

m.mad.traj <- ggplot(traj.sub, aes(x = lambda, y = exprs.diff.med)) + geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Pseudotime") + ylab("Median Log2 FC from Prog Cell") + 
  ggtitle(jmark, paste(trajname, "Genome-wide"))
print(m.mad.traj)



# H3K4me3 -----------------------------------------------------------------


jmark <- "H3K4me3"
trajname <- "myeloid"

print(paste(jmark, trajname))

imputed.dat <- t(tm.result.lst[[jmark]]$terms) %*% t(tm.result.lst[[jmark]]$topics)
dat.umap.long <- dat.umap.long.trajs[[jmark]]



cell.sd.df.long <- lapply(grep.strs, function(grep.str){
  return(GetCellSd(jscale * (imputed.dat + pseudo), grep.str, log2.scale = TRUE, fn = sd))
}) %>%
  bind_rows()

dat.umap.filt <- left_join(dat.umap.long, cell.sd.df.long)

print(range(dat.umap.filt$cell.sd))

m.chr <- ggplot(dat.umap.filt, aes(x = umap1, y = umap2, color = cell.sd)) + 
  geom_point(size = jsize.facet)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.sd, limits = lims.sd) + 
  ggtitle(jmark, paste0(deparse(substitute(sd)), " across chromosome"))
print(m.chr)

# what about genome wide
cell.sd.genomewide <- GetCellSd(imputed.dat, "", log2.scale=TRUE)
dat.umap.filt.gw <- left_join(dat.umap.long, cell.sd.genomewide)

m.gw <- ggplot(dat.umap.filt.gw, aes(x = umap1, y = umap2, color = cell.sd)) + 
  geom_point(size = gw.jsize.facet)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~label) + scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.sd, limits = lims.sd) + 
  ggtitle(jmark, paste0(deparse(substitute(sd)), " genome wide"))

print(m.gw)


# Highlight differences for two representative cells on a representative chromosome 
hsc.cell <- (trajs[[jmark]][[trajname]] %>% arrange(lambda) %>% dplyr::top_n(-1))$cell[[1]]
diff.cell <- (trajs[[jmark]][[trajname]] %>% arrange(lambda) %>% dplyr::top_n(1))$cell[[1]]

gstr <- paste0("chr15:")
jsub <- MatToLong(imputed.dat, gstr, cells.vec = c(hsc.cell, diff.cell))
m.spatial <- ggplot(jsub, aes(x = pos / 10^6, y = log2(exprs))) + 
  geom_line(alpha = jalpha) + 
  facet_wrap(~cell) + 
  ggtitle(paste(jmark, gstr)) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
m.spatial.merged <- ggplot(jsub, aes(x = pos / 10^6, y = log2(exprs), group = cell, color = cell)) + 
  geom_line(alpha = jalpha) + 
  ggtitle(paste(jmark, gstr)) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
m.spatial.log2fc <- ggplot(jsub %>% group_by(pos) %>% summarise(exprs = diff(log2(exprs))), aes(x = pos / 10^6, y = exprs)) + 
  geom_line(alpha = jalpha) + 
  ggtitle(paste(jmark, gstr)) + 
  xlab("MB") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ylab("log2 Fold Change")

print(m.spatial)
print(m.spatial.merged)
print(m.spatial.log2fc)

# spatial pattern?

jsub.hsc <- jsub %>% filter(cell == hsc.cell) %>% arrange(pos)
jsub.myeloid <- jsub %>% filter(cell == diff.cell) %>% arrange(pos)


# impute missing bins with minimum value

# jmain <- "hsc"
acf.out.hsc <- CalculateACF(jsub.hsc, jstep = jstep, jtype = jtype, jmain = paste("Progenitor Cell", gstr), show.plot = TRUE)
acf.out.myeloid <- CalculateACF(jsub.myeloid, jstep = jstep, jtype = jtype, jmain = paste("Differentiated Cell", gstr), show.plot = TRUE)

# do it genome wide
jsub.hsc.lst <- lapply(grep.strs, function(g) MatToLong(imputed.dat, g, cells.vec = c(hsc.cell)))
jsub.myeloid.lst <- lapply(grep.strs, function(g) MatToLong(imputed.dat, g, cells.vec = c(diff.cell)))

act.out.hsc.lst <- lapply(jsub.hsc.lst, function(jsub.hsc) CalculateACF(jsub.hsc, jstep = jstep, jtype = jtype, jmain = paste("Progenitor Cell", gstr), show.plot = FALSE))
act.out.myeloid.lst <- lapply(jsub.myeloid.lst, function(jsub.hsc) CalculateACF(jsub.hsc, jstep = jstep, jtype = jtype, jmain = paste("Differentiated Cell", gstr), show.plot = FALSE))

# average out the plots for different lags 
plot(act.out.hsc.lst[[1]]$lag.stepadj, act.out.hsc.lst[[1]]$acf)


# Plot the median log2 fold change relative to HSC cell: for one chromo

jsub.ref.merge <- lapply(grep.strs, function(gstr) GetDiffRelToCell(imputed.dat, gstr = gstr, trajs, trajname = trajname, dat.umap.long)) %>%
  bind_rows() 
m.mad <- ggplot(jsub.ref.merge, aes(x = umap1, y = umap2, color = exprs.diff.med)) + geom_point(size = jsize.facet) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.fc, limits = lims.fc) + facet_wrap(~label) + ggtitle(jmark, "Median log2 fold change with HSC")
print(m.mad)


# do genome-wide?
jsub.ref.merge.gw <- lapply(c(""), function(gstr) GetDiffRelToCell(imputed.dat, gstr = gstr, trajs, trajname = trajname, dat.umap.long)) %>%
  bind_rows() 
m.mad.gw <- ggplot(jsub.ref.merge.gw, aes(x = umap1, y = umap2, color = exprs.diff.med)) + geom_point(size = gw.jsize.facet) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.fc, limits = lims.fc) + facet_wrap(~label) + ggtitle(jmark, "Median log2 fold change with HSC")
print(m.mad.gw)

print(range(jsub.ref.merge.gw$exprs.diff.med))
# plot along pseudotime? 
traj.sub <- trajs[[jmark]][[trajname]]
# add exprs.diff.med
traj.sub <- left_join(traj.sub, jsub.ref.merge.gw %>% dplyr::select(cell, exprs.diff.med), by = c("cell"))

m.mad.traj <- ggplot(traj.sub, aes(x = lambda, y = exprs.diff.med)) + geom_point(alpha = 0.1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Pseudotime") + ylab("Median Log2 FC from Prog Cell") + 
  ggtitle(jmark, "Genome-wide")
print(m.mad.traj)



# do genome-wide?
jsub.ref.merge.gw <- lapply(c(""), function(gstr) GetDiffRelToCell(imputed.dat, gstr = gstr, trajs, trajname = trajname, dat.umap.long)) %>%
  bind_rows() 
m.mad.gw <- ggplot(jsub.ref.merge.gw, aes(x = umap1, y = umap2, color = exprs.diff.med)) + geom_point(size = gw.jsize.facet) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_gradient2(low = jcol[[1]], mid = jcol[[2]], high = jcol[[3]], midpoint = mdpt.fc, limits = lims.fc) + facet_wrap(~label) + ggtitle(jmark, "Median log2 fold change with HSC")
print(m.mad.gw)

print(range(jsub.ref.merge.gw$exprs.diff.med))
# plot along pseudotime? 
traj.sub <- trajs[[jmark]][[trajname]]
# add exprs.diff.med
traj.sub <- left_join(traj.sub, jsub.ref.merge.gw %>% dplyr::select(cell, exprs.diff.med), by = c("cell"))

m.mad.traj <- ggplot(traj.sub, aes(x = lambda, y = exprs.diff.med)) + geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("Pseudotime") + ylab("Median Log2 FC from Prog Cell") + 
  ggtitle(jmark, paste(trajname, "Genome-wide"))
print(m.mad.traj)



if (do.plots){
  dev.off()
}

print(Sys.time() - tstart)