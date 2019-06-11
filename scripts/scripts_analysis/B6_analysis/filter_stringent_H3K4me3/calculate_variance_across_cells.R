# Jake Yeung
# Date of Creation: 2019-06-04
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/filter_stringent_H3K4me3/calculate_variance_across_cells.R
# Load LDA, calculate variance 

rm(list=ls())
library(dplyr)
library(ggplot2)
library(topicmodels)
library(umap)
library(Matrix)

source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/PlotFunctions.R")


# Functions ---------------------------------------------------------------

FitGetPval <- function(jsub, jform){
  jfit <- lm(formula = jform, data = jsub)
  pval <- summary(jfit)$coefficients[2, "Pr(>|t|)"]
  slope.val <- summary(jfit)$coefficients[2, "Estimate"]
  slope.se <- summary(jfit)$coefficients[2, "Std. Error"]
  return(data.frame(pval = pval, slope.val = slope.val, slope.se = slope.se))
}

SumSqrDev <- function(x){
  return( sum( (x - mean(x)) ^ 2 ))
}

# Load data ---------------------------------------------------------------



inf <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.H3K4me3.stringent_filter.RData"
load(inf, v=T)
tm.result.stringent <- posterior(out.objs$out.lda)
dat.umap.long.stringent <- dat.umap.long

inf.cellsums <- "/Users/yeung/data/scchic/robjs/B6_objs/cell_sums_from_LDA_input.RData"
load(inf.cellsums, v=T)

inf.objs <- "/Users/yeung/data/scchic/robjs/B6_objs_stringent/traj_objs_H3K4me3.stringent.Rdata"
load(inf.objs, v=T)
dat.umap.long.trajs.stringent <- dat.umap.long.trajs
trajs.stringent <- trajs
trajs.objs.stringent <- trajs.objs
# trajs.long.stringent <- trajs.long

jfac <- 10^6
jpseudo <- 0
dat.mat <-  t(tm.result.stringent$terms) %*% t(tm.result.stringent$topics)
# log2 transform
dat.mat <- log2(dat.mat * jfac + jpseudo)

cells.sd <- GetCellSd(dat.mat, "", log2.scale = FALSE, fn = var) 

# merge with umap and plot

dat.merged.stringent <- left_join(dat.umap.long.stringent, cells.sd) %>%
  dplyr::rename(cell.var = cell.sd)

PlotXYWithColor(dat.merged.stringent, xvar = "umap1", yvar = "umap2", cname = "cell.var", jsize = 5)


# Do intra and inter chromosomal variance  --------------------------------



# Add cell counts along trajectory  ---------------------------------------

dat.sum.merged <- left_join(dat.umap.long.stringent, count.sum)
# add trajs info
jtrajs <- c("granu", "lymphoid", "eryth")
trajs.long.stringent <- lapply(jtrajs, function(jtraj){
  trajs.tmp <- lapply(trajs.stringent, function(x) x[[jtraj]]) %>%
    bind_rows() %>%
    rowwise() %>%
    mutate(mark = strsplit(cell, "-")[[1]][[4]]) %>%
    rowwise() %>%
    mutate(lambda.bin = floor(lambda * 10) / 10) %>%
    mutate(traj = jtraj)
}) %>%
  bind_rows()

# add cell counts to trajs.long and plot

trajs.long.merged <- left_join(trajs.long.stringent, count.sum)

trajname <- "granu"

ggplot(trajs.long.merged, aes(x = lambda, y = log10(cellsum), color = traj)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Regenerate FIgure 4  ----------------------------------------------------


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

# add stringent
tm.result.lst$H3K4me3 <- tm.result.stringent

# Load trajectories -------------------------------------------------------

inf.traj <- paste0("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata")
load(inf.traj, v=T)
dat.umap.long.trajs$H3K4me3 <- dat.umap.long.trajs.stringent$H3K4me3
trajs$H3K4me3 <- trajs.stringent$H3K4me3
trajs.objs$H3K4me3 <- trajs.objs.stringent$H3K4me3

# Calculate variance per cell  --------------------------------------------

cells.sd <- lapply(jmarks, function(jmark){
  dat.mat <-  t(tm.result.lst[[jmark]]$terms) %*% t(tm.result.lst[[jmark]]$topics)
  # log2 transform
  dat.mat <- log2(dat.mat * jfac + jpseudo)
  cells.sd <- GetCellSd(dat.mat, "", log2.scale = FALSE) %>%
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

# pdf(paste0("/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory/variance_over_pseudotime_fewer_trajs.", Sys.Date(), ".pdf"), useDingbats = FALSE)
print(m.facet)
print(m.nofacet)
# dev.off()

# Break down variance by chromosomes  -------------------------------------


# Break down variance by chromosomes  -------------------------------------



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
    # jtest <- apply(dat.mat, MARGIN = 2, FUN = SumSqrDev)
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

cells.var.long <- cells.var.merged %>%
  dplyr::select(-cell.var.across.expect) %>%
  tidyr::gather(key = varname, value = varval, -c(cell, label, mark, jdiff, within.frac))
cells.var.merged2 <- left_join(cells.var.long, trajs.long %>% dplyr::select(mark, cell, lambda, lambda.bin, traj))
cells.var.merged2$mark <- factor(cells.var.merged2$mark, levels = c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"))

# normalization factor
chromo.constant.sum <- chromo.constant %>%
  summarise(nbins = sum(nbins))

jtrajname <- "granu"

# save outputs
save(cells.var.merged2, cells.chromo.mean.wide.mat, dat.mat.global.mean, chromo.constant, cells.var.chromo.within, file = "/Users/yeung/data/scchic/robjs/B6_objs_stringent/variance_decomposed_within_across_stringent.RData")

pdf(paste0("/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory/decompose_variance.", Sys.Date(), ".stringent.pdf"), useDingbats = FALSE)
for (jtrajname in jtrajs){
  m.var <- ggplot(cells.var.merged2 %>% filter(traj == jtrajname), aes(x = lambda, y = varval / chromo.constant.sum$nbins, color = varname, group = varname)) + geom_point(alpha = 0.2) + facet_wrap(~mark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_smooth(method = "lm", se = FALSE) + 
    xlab("Pseudotime") + ylab("Variance (log2 signal ^ 2)") + 
    scale_color_discrete(name = "", labels = c("VarGlobal", "VarAcross", "VarWithin")) + 
    ggtitle(jtrajname)
  print(m.var)
  
  ssdiff <- sweep(x = cells.chromo.mean.wide.mat, MARGIN = 2, STATS = dat.mat.global.mean, FUN = "-") ^ 2
  ssdiff <- sweep(x = ssdiff, MARGIN = 1, STATS = chromo.constant$nbins, FUN = "*")
  chromo.mean.dat <- data.frame(chromo.var = unlist(data.frame(ssdiff)), 
                                chromo = rep(rownames(ssdiff), ncol(ssdiff)), 
                                cell = rep(colnames(ssdiff), each = nrow(ssdiff)), stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(mark = strsplit(cell, "-")[[1]][[4]])
  m.withinbychromo <- ggplot(cells.var.chromo.within, aes(x = log10(cell.var.within), fill = mark)) + geom_density(alpha = 0.25) + facet_wrap(~label)
  m.acrossvar <- ggplot(chromo.mean.dat, aes(x = log10(chromo.var), fill = mark)) + geom_density(alpha = 0.25) + facet_wrap(~chromo) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m.withinbychromo)
  print(m.acrossvar)
}
dev.off()




# pdf(pdfout, useDingbats = FALSE)

pdf(paste0("/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory/variance_over_pseudotime_fewer_trajs.withinvar.", Sys.Date(), ".stringent.pdf"), useDingbats=FALSE)

# variance over pseudotime 
colhash <- GetTrajColors(as.hash = TRUE, add.mega = FALSE)
jsub <- cells.var.merged2 %>% filter(!is.na(traj))
jsub$jcol <- sapply(jsub$traj, function(x) colhash[[x]])
m.var1 <- ggplot(jsub %>% filter(varname == "cell.var.within.sum"), aes(x = lambda, y = varval / chromo.constant.sum$nbins, color = jcol)) + geom_point(alpha = 0.2) + facet_wrap(~mark, nrow = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_continuous(breaks = c(0, 1)) + 
  xlab("Pseudotime") + 
  ylab(expression('Variance ' * group("[", log[2] * ("signal")^{2}, "]"))) + 
  scale_color_identity()
print(m.var1)
# dev.off()

m.var.smooth <- ggplot(jsub %>% filter(varname == "cell.var.within.sum"), aes(x = lambda, y = varval / chromo.constant.sum$nbins, color = jcol)) + geom_point(alpha = 0.2) + facet_wrap(~mark, nrow = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_smooth(method = "lm", se = FALSE) + 
  scale_x_continuous(breaks = c(0, 1)) + 
  xlab("Pseudotime") +
  ylab(expression('Variance ' * group("[", log[2] * ("signal")^{2}, "]"))) + 
  scale_color_identity()
print(m.var.smooth)

# umap
dat.merged <- left_join(bind_rows(dat.umap.long.trajs), cells.var.merged2 %>% filter(varname == "cell.var.within.sum") %>% dplyr::select(cell, varval)) %>%
  mutate(cell.within.var = varval / chromo.constant.sum$nbins)

# library("colorspace")
# pal <- colorspace::choose_palette()

jtitle <- ""
jcname <- "cell.within.var"
xvar <- "umap1"
yvar <- "umap2"
# steps <- c("blue4", "cyan", "white", "yellow", "red4")
steps <- c("gray85", "gray50", "gray40",  "cyan")
pal <- color.palette(steps, c(6, 2, 40), space="rgb")
jsize <- 4

for (jmark in jmarks){
  # m1 <- PlotXYWithColor(dat.merged %>% filter(mark == jmark), xvar = "umap1", yvar = "umap2", cname = "cell.within.var", jcol.low = "gray90", jcol.mid = "gray50", jcol = "cyan", jsize = 6)
  # print(m1)
  dat.merged$neg.cell.within.var <- -1 * dat.merged$cell.within.var
  dat.merged <- RankOrder(dat.merged, cname = "neg.cell.within.var", out.cname = "orderrank")
  m1 <- ggplot(dat.merged %>% filter(mark == jmark), aes_string(x = xvar, y = yvar, col = jcname, order = "orderrank")) +
    ggrastr::geom_point_rast(size = jsize) +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
                       axis.ticks=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank(),
                       panel.border=element_blank())  +
    xlab("") + ylab("") + ggtitle(jmark) +
    # scale_color_gradientn(colours = RColorBrewer::brewer.pal(20, "YlGnBu"))
    scale_color_gradientn(colours = pal(50))
  # scale_color_gradient2(low = "gray85", high = "cyan", mid = "gray50", breaks = c(0.5, 2, 4), midpoint = 2)
  print(m1)
}

for (jtrajname in jtrajs){
  m.var <- ggplot(cells.var.merged2 %>% filter(traj == jtrajname), aes(x = lambda, y = varval / chromo.constant.sum$nbins, color = varname, group = varname)) + geom_point(alpha = 0.2) + facet_wrap(~mark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_smooth(method = "lm", se = FALSE, size = 2) + 
    xlab("Pseudotime") + ylab("Variance (log2 signal ^ 2)") + 
    scale_color_discrete(name = "", labels = c("VarGlobal", "VarAcross", "VarWithin")) + 
    ggtitle(jtrajname)
  print(m.var)
  
  ssdiff <- sweep(x = cells.chromo.mean.wide.mat, MARGIN = 2, STATS = dat.mat.global.mean, FUN = "-") ^ 2
  ssdiff <- sweep(x = ssdiff, MARGIN = 1, STATS = chromo.constant$nbins, FUN = "*")
  chromo.mean.dat <- data.frame(chromo.var = unlist(data.frame(ssdiff)), 
                                chromo = rep(rownames(ssdiff), ncol(ssdiff)), 
                                cell = rep(colnames(ssdiff), each = nrow(ssdiff)), stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(mark = strsplit(cell, "-")[[1]][[4]])
  m.withinbychromo <- ggplot(cells.var.chromo.within, aes(x = log10(cell.var.within), fill = mark)) + geom_density(alpha = 0.25) + facet_wrap(~label)
  m.acrossvar <- ggplot(chromo.mean.dat, aes(x = log10(chromo.var), fill = mark)) + geom_density(alpha = 0.25) + facet_wrap(~chromo) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m.withinbychromo)
  print(m.acrossvar)
}
dev.off()


# Do UMIs over differentiation  -------------------------------------------



count.sum <- count.sum %>%
  dplyr::rename(cell.sum = cellsum) %>%
  mutate(cell.sum.log10 = log10(cell.sum))

cell.sums <- left_join(count.sum, dat.umap.long.trajs %>% bind_rows() %>% dplyr::select(cell, umap1, umap2, mark))

traj.info.all <- lapply(jmarks, function(jmark){
  traj.info <- lapply(jtrajs, function(jtraj){
    return(trajs[[jmark]][[jtraj]] %>% mutate(traj = jtraj, mark = jmark) %>% dplyr::select(c(-umap1, -umap2)))
  }) %>%
    bind_rows()
}) %>% bind_rows()

cell.sums.merged.all <- left_join(cell.sums, traj.info.all)
cell.sums.merged.all$jcol <- sapply(cell.sums.merged.all$traj, function(x) colhash[[x]])
cell.sums.merged.all$mark <- factor(as.character(cell.sums.merged.all$mark), levels = jmarks)

cell.sums$neg.sum <- -1 * cell.sums$cell.sum
cell.sums<- RankOrder(cell.sums, cname = "neg.sum", out.cname = "orderrank")
cell.sums<- RankOrder(cell.sums, cname = "cell.sum", out.cname = "orderrank2")

m.all.linear <- ggplot(cell.sums.merged.all %>% filter(!is.na(traj)), aes(x = lambda, y = cell.sum, color = jcol)) + geom_point(alpha = 0.3) + facet_wrap(~mark, nrow = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_x_continuous(breaks = c(0, 1), labels = c("", "")) + 
  ylab("cuts/cell") + 
  xlab("") + 
  scale_color_identity()

m.all <- ggplot(cell.sums.merged.all %>% filter(!is.na(traj)), aes(x = lambda, y = cell.sum, color = jcol)) + geom_point(alpha = 0.3) + facet_wrap(~mark, nrow = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_y_log10()  + 
  scale_x_continuous(breaks = c(0, 1), labels = c("", "")) + 
  ylab("cuts/cell") + 
  xlab("") + 
  scale_color_identity()

m.all.log <- ggplot(cell.sums.merged.all %>% filter(!is.na(traj)), aes(x = lambda, y = cell.sum.log10, color = jcol)) + geom_point(alpha = 0.3) + facet_wrap(~mark, nrow = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_continuous(breaks = c(0, 1), labels = c("", "")) + 
  ylab(expression(log[10]~"(cuts/cell)")) + 
  xlab("") + 
  scale_color_identity()

m.all.log.smooth <- ggplot(cell.sums.merged.all %>% filter(!is.na(traj)), aes(x = lambda, y = cell.sum.log10, color = jcol, group = traj)) + geom_point(alpha = 0.3) + facet_wrap(~mark, nrow = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_continuous(breaks = c(0, 1), labels = c("", "")) + 
  ylab(expression(log[10]~"(cuts/cell)")) + 
  xlab("") + 
  geom_smooth(method = "lm", se = FALSE) + 
  scale_color_identity()

steps <- c("gray95", "gray70", "gray60",  "red")
pal.sums <- color.palette(steps, c(10, 5, 10), space="rgb")

# pdf(paste0("/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory/cell_counts_over_umap.", Sys.Date(), ".stringent.pdf"), useDingbats=FALSE)
pdf(paste0("/Users/yeung/data/scchic/pdfs/B6_figures/stringent_pdfs/variance_over_trajectory_stringent.", Sys.Date(), ".stringent.pdf"), useDingbats=FALSE)
# plot cell sums in linear
jcname <- "cell.sum.log10"
for (jmark in jmarks){
  if (jmark == "H3K9me3"){
    ordervar <- "orderrank2"
  } else {
    ordervar <- "orderrank"
  }
  m1 <- ggplot(cell.sums %>% filter(mark == jmark), aes_string(x = xvar, y = yvar, col = jcname, order = ordervar)) +
    ggrastr::geom_point_rast(size = jsize) +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
                       axis.ticks=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank(),
                       panel.border=element_blank())  +
    xlab("") + ylab("") + ggtitle(jmark) +
    # scale_color_gradientn(colours = RColorBrewer::brewer.pal(9, "Reds"))
    scale_color_gradientn(colours = pal.sums(50))
  # scale_color_gradient2(low = "gray85", high = "cyan", mid = "gray50", breaks = c(0.5, 2, 4), midpoint = 2)
  print(m1)
  
  # m <- PlotXYWithColor(cell.sums %>% filter(mark == jmark), xvar = "umap1", yvar = "umap2", cname = "cell.sum", jcol = scales::muted("darkred"), jtitle = jmark, jsize = jsize)
  # print(m)
}
# plot cell sums in log 
for (jmark in jmarks){
  m.log <- PlotXYWithColor(cell.sums %>% filter(mark == jmark), xvar = "umap1", yvar = "umap2", cname = "cell.sum.log10", jcol = scales::muted("darkred"), jtitle = jmark, jsize = jsize)
  print(m.log)
}
print(m.all.linear)
print(m.all)
print(m.all.log)
print(m.all.log.smooth)

for (jmark in jmarks){
  traj.info <- lapply(jtrajs, function(jtraj){
    return(trajs[[jmark]][[jtraj]] %>% mutate(traj = jtraj, mark = jmark) %>% dplyr::select(c(-umap1, -umap2)))
  }) %>%
    bind_rows()
  cell.sums.merged <- left_join(cell.sums, traj.info)
  cell.sums.merged$jcol <- sapply(cell.sums.merged$traj, function(x) ifelse(!is.null(colhash[[x]]), colhash[[x]], NA))
  m.traj <- ggplot(cell.sums.merged %>% filter(!is.na(traj)), aes(x = lambda, y = cell.sum.log10, color = jcol)) + facet_wrap(~traj) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_identity() + 
    ggtitle(jmark)
  m.traj.lin <- ggplot(cell.sums.merged %>% filter(!is.na(traj)), aes(x = lambda, y = cell.sum, color = jcol)) + facet_wrap(~traj) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_identity() + 
    ggtitle(jmark)
  print(m.traj)
  print(m.traj.lin)
}
# put marks in one graphs

# jtmp <- cell.sums.merged.all %>% filter(!is.na(traj)) %>% filter(mark == "H3K4me3" & traj == "granu")
# jfit.tmp <- lm(formula = as.formula("cell.sum.log10 ~ lambda"), jtmp)


# fit the intrachromosomal variance over pseudotime

# ggplot(jfits, aes(x = traj, y = -log10(pval))) + facet_wrap(~mark) + geom_bar(stat = "identity") + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jfits <- cell.sums.merged.all %>%
  filter(!is.na(traj)) %>%
  group_by(mark, traj) %>%
  mutate(cell.sum.log2 = log2(cell.sum)) %>%
  do(FitGetPval(., jform = as.formula("cell.sum.log2 ~ lambda"))) %>% 
  rowwise() %>%
  mutate(jcol = colhash[[traj]]) %>%
  arrange(desc(pval))

jfits.merged <- cell.sums.merged.all %>%
  filter(!is.na(traj)) %>%
  group_by(mark) %>%
  mutate(cell.sum.log2 = log2(cell.sum)) %>%
  do(FitGetPval(., jform = as.formula("cell.sum.log2 ~ lambda"))) %>%
  rowwise() %>%
  mutate(jcol = colhash[[traj]]) %>%
  arrange(desc(pval))

jfits.var <- cells.var.merged2 %>% filter(!is.na(traj)) %>%
  group_by(mark, traj) %>%
  filter(varname == "cell.var.within.sum") %>%
  mutate(varval = varval / chromo.constant.sum$nbins) %>%
  do(FitGetPval(., jform = as.formula("varval ~ lambda"))) %>%
  rowwise() %>%
  mutate(jcol = colhash[[traj]])

jfits.merged.var <- cells.var.merged2 %>% filter(!is.na(traj)) %>%
  group_by(mark) %>%
  filter(varname == "cell.var.within.sum") %>%
  mutate(varval = varval / chromo.constant.sum$nbins) %>%
  do(FitGetPval(., jform = as.formula("varval ~ lambda"))) %>%
  rowwise() %>%
  mutate(jcol = colhash[[traj]])

m.reg.sum <- ggplot(jfits, aes(x = traj, y = slope.val, ymin = slope.val - slope.se, ymax = slope.val + slope.se, fill = jcol)) + 
  facet_wrap(~mark, nrow = 1) + 
  geom_bar(stat = "identity") + geom_errorbar(width = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_identity() + 
  xlab("") + ylab("Log2 FC in # cuts")

m.reg.var <- ggplot(jfits.var, aes(x = traj, y = slope.val, ymin = slope.val - slope.se, ymax = slope.val + slope.se, fill = jcol)) + 
  facet_wrap(~mark, nrow = 1) + 
  geom_bar(stat = "identity") + geom_errorbar(width = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_identity() + 
  xlab("") + ylab("Change in Variance")

print(m.reg.sum)
print(m.reg.var)
dev.off()


# Do linear model to get p-value of the slope  ----------------------------

# jsub <- cell.sums.merged.all %>% filter(!is.na(traj)) %>% filter(mark == "H3K9me3" & traj == "eryth")
# 
# 
# 
# jtmp <- cell.sums.merged.all %>% filter(!is.na(traj)) %>% filter(mark == "H3K4me3" & traj == "granu")
# jfit.tmp <- lm(formula = as.formula("cell.sum.log10 ~ lambda"), jtmp)
# 
# (jpval <- FitGetPval(cell.sums.merged.all %>% filter(!is.na(traj)) %>% filter(mark == "H3K4me3" & traj == "granu"), jform = as.formula("cell.sum.log10 ~ lambda")))
# (jpval <- FitGetPval(cell.sums.merged.all %>% filter(!is.na(traj)) %>% filter(mark == "H3K4me1" & traj == "lymphoid")))
# (jpval <- FitGetPval(cell.sums.merged.all %>% filter(!is.na(traj)) %>% filter(mark == "H3K4me1" & traj == "eryth")))

# show slope value with error bar 




# Save tables  ------------------------------------------------------------


save(jfits, jfits.merged, cell.sums.merged.all, file = "/Users/yeung/data/scchic/robjs/fits_cell_sum_along_traj/fits_cell_size_along_traj.RData")
save(jfits.var, jfits.merged.var, cells.var.merged2, file = "/Users/yeung/data/scchic/robjs/fits_cell_sum_along_traj/variance_along_traj.RData")

# write tables to dropbox
inf1 <- "/Users/yeung/Dropbox/scCHiC_figs/tables/fits_counts_along_traj.txt"
inf2 <- "/Users/yeung/Dropbox/scCHiC_figs/tables/fits_merged_counts_along_traj.txt"
inf3 <- "/Users/yeung/Dropbox/scCHiC_figs/tables/fits_var_along_traj.txt"
inf4 <- "/Users/yeung/Dropbox/scCHiC_figs/tables/fits_merged_var_along_traj.txt"

write.table(jfits, file = inf1, quote = FALSE, sep = "\t", row.names = FALSE)
write.table(jfits.merged, file = inf2, quote = FALSE, sep = "\t", row.names = FALSE)
write.table(jfits.var, file = inf3, quote = FALSE, sep = "\t", row.names = FALSE)
write.table(jfits.merged.var, file = inf4, quote = FALSE, sep = "\t", row.names = FALSE)


