# Jake Yeung
# Date of Creation: 2019-05-14
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/8-variance_over_pseudotime.R
# Variance over pseudotime break down by chromosome and also across chromosome 

rm(list=ls())

library(data.table)
library(dplyr)
library(ggplot2)
library(JFuncs)
library(tidyr)
library(GGally)
library(ggrastr)
library(gridExtra)

library(grDevices)

library(Matrix)

source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/PlotFunctions.R")

# Get constants -----------------------------------------------------------

make.plots <- TRUE
colhash <- GetTrajColors(as.hash = TRUE, add.mega = FALSE)
jtrajs <- c("granu", "lymphoid", "eryth")

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
colvec <- c("gray85", "gray50", "blue")  
names(jmarks) <- jmarks

nsecs <- 5
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

# do.plots <- TRUE

jstep <- 20000
# jtype <- "correlation"

pos.max <- 50 * 10^6
jstep <- 20000
lagmax <- pos.max / jstep

ci.interval <- 1.96  # corresponds to zscore 1.96 = alpha = 0.05 = 95% confidence interval

jsize <- 2
jthres <- 0.05
jpseudo <- 0
jfac <- 10^6

jtraj <- "granu"

# Load data  --------------------------------------------------------------

indir <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmark <- jmarks[["H3K4me1"]]
infs <- lapply(jmarks, function(jmark) list.files(indir, pattern = paste0(jmark, ".RData"), full.names = TRUE))
tm.result.lst <- lapply(infs, LoadGetTmResult)


# Load trajectories -------------------------------------------------------

inf.traj <- paste0("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata")
load(inf.traj, v=T)


# Calculate variance per cell  --------------------------------------------

cells.sd <- lapply(jmarks, function(jmark){
  dat.mat <-  t(tm.result.lst[[jmark]]$terms) %*% t(tm.result.lst[[jmark]]$topics)
  # log2 transform
  dat.mat <- log2(dat.mat * jfac + jpseudo)
  cells.sd <- GetCellSd(dat.mat, "", log2.scale = FALSE, fn = sd) %>%
    mutate(mark = jmark)
  return(cells.sd)
})
cells.sd <- cells.sd %>%
  bind_rows()

# add lambda
trajs.long <- lapply(jtrajs, function(jtraj){
  trajs.tmp <- lapply(trajs, function(x) x[[jtraj]]) %>%
    bind_rows() %>%
    rowwise() %>%
    mutate(mark = strsplit(cell, "-")[[1]][[4]]) %>%
    rowwise() %>%
    mutate(lambda.bin = floor(lambda * 10) / 10) %>%
    mutate(traj = jtraj)
}) %>%
  bind_rows()

# add info
cells.sd.merge <- left_join(cells.sd, trajs.long %>% dplyr::select(mark, cell, lambda, lambda.bin, traj))
cells.sd.merge$mark <- factor(cells.sd.merge$mark, levels = c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"))


jsub <- cells.sd.merge %>% filter(!is.na(traj)) %>% rowwise() %>% mutate(jcol = colhash[[traj]])
jsub.trajfilt <- cells.sd.merge %>% filter(traj == jtraj) %>% rowwise() %>% mutate(jcol = colhash[[jtraj]])
m.facet <- ggplot(jsub, aes(x = lambda, y = cell.sd, color = jcol, group = traj)) + 
  # facet_wrap(~mark, nrow = 1) + 
  facet_grid(traj~mark) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text.x = element_blank(), axis.ticks.x = element_blank()) +  
  scale_color_identity() + 
  xlab("Pseudotime") + ylab("Genome-wide SD") 

m.nofacet <- ggplot(jsub, aes(x = lambda, y = cell.sd, color = jcol, group = traj)) + 
  facet_wrap(~mark, nrow = 1) + 
  geom_point(alpha = 0.3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                     axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom") +  
  scale_color_identity() + 
  xlab("Pseudotime") + ylab("Genome-wide SD") 

pdf(paste0("/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory/variance_over_pseudotime_fewer_trajs.", Sys.Date(), ".pdf"), useDingbats = FALSE)
  print(m.facet)
  print(m.nofacet)
dev.off()

# Break down variance by chromosomes  -------------------------------------

SumSqrDev <- function(x){
  return( sum( (x - mean(x)) ^ 2 ))
}

cells.var <- lapply(jmarks, function(jmark){
  dat.mat <-  t(tm.result.lst[[jmark]]$terms) %*% t(tm.result.lst[[jmark]]$topics)
  # log2 transform
  dat.mat <- log2(dat.mat * jfac + jpseudo)
  # jtest <- apply(dat.mat, MARGIN = 2, FUN = SumSqrDev)
  cells.sd <- GetCellSd(dat.mat, "", log2.scale = FALSE, fn = SumSqrDev) %>%
    mutate(mark = jmark)
  return(cells.sd)
}) %>%
  bind_rows() %>%
  dplyr::rename(cell.var = cell.sd)

# do variance across chromosomes
jchromos <- c(paste("chr", seq(19), sep = ""), "chrX", "chrY")
jchromos.grep <- paste(jchromos, ":", sep = "")
jmark <- "H3K4me1"

# calculate variance within chromosomes
cells.var.chromo.within <- lapply(jmarks, function(jmark){
  cells.var.chromo <- lapply(jchromos.grep, function(jstr){
    dat.mat <-  t(tm.result.lst[[jmark]]$terms[, grepl(jstr, colnames(tm.result.lst[[jmark]]$terms))]) %*% t(tm.result.lst[[jmark]]$topics)
    # log2 transform
    dat.mat <- log2(dat.mat * jfac + jpseudo)
    cells.sd <- GetCellSd(dat.mat, grep.str = jstr, log2.scale = FALSE, fn = SumSqrDev) %>%
      mutate(mark = jmark)
    jtest <- apply(dat.mat, MARGIN = 2, FUN = SumSqrDev)
    return(cells.sd)
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  dplyr::rename(cell.var.within = cell.sd)

# summarize across chromosomes
cells.var.chromo.within.sum <- cells.var.chromo.within %>%
  group_by(cell) %>%
  summarise(cell.var.within.sum = sum(cell.var.within))

# calculate variance across chromosomes
cells.chromo.mean <- lapply(jmarks, function(jmark){
  cells.chromo.mean.tmp <- lapply(jchromos.grep, function(jstr){
    dat.mat <-  t(tm.result.lst[[jmark]]$terms[, grepl(jstr, colnames(tm.result.lst[[jmark]]$terms))]) %*% t(tm.result.lst[[jmark]]$topics)
    # log2 transform
    dat.mat <- log2(dat.mat * jfac + jpseudo)
    # get average across chromosomes
    chromo.avg <- colMeans(dat.mat)
    # chromo.avg <- colSums(dat.mat)
    chromo.avg.dat <- data.frame(cell = names(chromo.avg), chromo.mean = chromo.avg, label = jstr, nbins = nrow(dat.mat))
    return(chromo.avg.dat)
  }) %>%
    bind_rows()
}) %>%
  bind_rows()

# calculate global mean
dat.mat.global.mean <- lapply(jmarks, function(jmark){
  dat.mat.global <- log2(t(tm.result.lst[[jmark]]$terms) %*% t(tm.result.lst[[jmark]]$topics) * jfac + jpseudo)
  out <- colMeans(dat.mat.global)
  return(out)
}) %>%
  unlist()
# make sure cell names align with cells.chromo.mean.wide.mat
names(dat.mat.global.mean) <- sapply(names(dat.mat.global.mean), function(x) strsplit(x, "\\.")[[1]][[2]])
dat.mat.global.mean <- dat.mat.global.mean[order(names(dat.mat.global.mean))]

chromo.constant <- cells.chromo.mean %>%
  group_by(label) %>%
  summarise(nbins = unique(nbins))
  
cells.chromo.mean.wide <- tidyr::spread(cells.chromo.mean, key = cell, value = chromo.mean)
cells.chromo.mean.wide.mat <- as.matrix(cells.chromo.mean.wide %>% dplyr::select(-label, -nbins))
rownames(cells.chromo.mean.wide.mat) <- cells.chromo.mean.wide$label
# sum across chromosomes
# cells.var.chromo.across <- colSums(sweep(x = cells.chromo.mean.wide.mat, MARGIN = 2, STATS = dat.mat.global.mean, FUN = "-") ^ 2)

# multply constant for each chromosome
assertthat::assert_that(all(names(dat.mat.global.mean) == colnames(cells.chromo.mean.wide.mat)))
cells.var.chromo.across <- sweep(x = cells.chromo.mean.wide.mat, MARGIN = 2, STATS = dat.mat.global.mean, FUN = "-") ^ 2
assertthat::assert_that(all(chromo.constant$label == rownames(cells.chromo.mean.wide.mat)))
cells.var.chromo.across <- sweep(x = cells.var.chromo.across, MARGIN = 1, STATS = chromo.constant$nbins, FUN = "*")
cells.var.chromo.across <- colSums(cells.var.chromo.across)

cells.var.chromo.across <- data.frame(cell = names(cells.var.chromo.across), cell.var.across = cells.var.chromo.across)

# check the addition


# merge cells.chromo.mean 
cells.var.merged <- left_join(cells.var.chromo.within.sum, cells.var.chromo.across)
cells.var.merged <- left_join(cells.var.merged, cells.var)

cells.var.merged$cell.var.across.expect <- cells.var.merged$cell.var - cells.var.merged$cell.var.within.sum
cells.var.merged$jdiff <- cells.var.merged$cell.var.across.expect - cells.var.merged$cell.var.across
cells.var.merged$within.frac <- cells.var.merged$cell.var.within.sum / cells.var.merged$cell.var
plot(density(cells.var.merged$jdiff))



# Break down variance by cells --------------------------------------------

plot(density(cells.var.merged$cell.var.within.sum / cells.var.merged$cell.var))

ggplot(cells.var.merged, aes(x = within.frac)) + geom_histogram() + facet_wrap(~mark)


# Add trajectory info and plot --------------------------------------------


