# Jake Yeung
# Date of Creation: 2019-05-03
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/variance_over_pseudotime_heatmaps_zoomed_in.R
# Zoom in on genomic regions 


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

library(GGally)

source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/Aux.R")

FitSlope <- function(dat.sub){
  # fit exprs to trajectory to find slope
  fit <- lm(formula = exprs ~ lambda.bin, data = dat.sub)
  slope <- fit$coefficients[['lambda.bin']]
  int <- fit$coefficients[['(Intercept)']]
  pval <- summary(fit)$coefficients["lambda.bin", "Pr(>|t|)"]
  return(data.frame(slope = slope, int = int, pval = pval))
}



# Make as function --------------------------------------------------------


jtraj <- "granu"
BinTrajectory <- function(trajs.spring.lst, jtraj, nearest = 0.1){
  round.int <- 1 / nearest
  trajs.sum <- lapply(trajs.spring, function(x) x[[jtraj]] %>% mutate(traj = jtraj)) %>%
    bind_rows() %>%
    rowwise() %>%
    mutate(mark = strsplit(cell, "_")[[1]][[2]]) %>%
    left_join(., mat.sub.merge) %>%
    rowwise() %>%
    mutate(lambda.bin = floor(lambda * round.int) / round.int) %>%
    group_by(traj, lambda.bin, mark, coord, pos) %>%
    summarise(exprs = mean(exprs)) %>%
    return(trajs.sum)
}



# Load data  --------------------------------------------------------------

# inf.trajs <- "/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-09.RData"
inf.trajs <- "/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-11.RData"
inf.trajs <- "/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-15.RData"
assertthat::assert_that(file.exists(inf.trajs))
load(inf.trajs, v=T)

head(trajs.spring[[4]][[1]])

inf.dat <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"
assertthat::assert_that(file.exists(inf.dat))
load(inf.dat, v=T)


# Get constants -----------------------------------------------------------

make.plots <- TRUE
make.plots <- FALSE

colvec <- c("gray85", "gray50", "blue")  
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
names(jmarks) <- jmarks
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#C0C0C0")

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


jstr <- "chr15:"

jthres <- 0.05

# Explore the data  -------------------------------------------------------

jmark <- "H3K4me1"
jmark <- "H3K27me3"
jmark <- "H3K9me3"

jtraj <- "granu"

# Prepare dat -------------------------------------------------------------

jpseudo <- 0
jfac <- 10^6


# imputed.dat <- log2(t(tm.result.lst[[jmark]]$terms) %*% t(tm.result.lst[[jmark]]$topics) + jpseudo) * jfac
# mat.sub <- MatToLong(imputed.dat, gstr = "chr15:", cells.vec = NULL) %>% dplyr::select(-start, -end)

mat.sub.merge <- lapply(jmarks, function(jmark) GetMatSub(tm.result.lst, jmark, jstr, jpseudo, jfac) %>% mutate(mark = jmark)) %>% 
  bind_rows()

head(trajs.spring[[1]][[1]])

jtraj <- "granu"
trajs.long <- lapply(trajs.spring, function(x) x[[jtraj]]) %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(mark = strsplit(cell, "_")[[1]][[2]]) %>%
  left_join(., mat.sub.merge) %>%
  rowwise() %>%
  mutate(lambda.bin = floor(lambda * 10) / 10)
trajs.sum <- trajs.long %>%
  group_by(lambda.bin, mark, coord, pos) %>%
  summarise(exprs = mean(exprs))

# Do heatmap of granulocytes for the 4 marks ------------------------------


jsub.fits <- trajs.sum %>%
  group_by(pos, mark) %>%
  do(FitSlope(.))

jsub.wide <- subset(jsub.fits, select = c(pos, mark, slope)) %>% spread(key = mark, value = slope) %>% ungroup()

chunksize <- 2e7  # 20 MB
trajs.sum$pos.round <- as.integer(floor(trajs.sum$pos / 2e7) * 2e7)
trajs.sum$pos.round2 <- as.integer(floor(trajs.sum$pos / 5e6) * 5e6)
trajs.sum$mark <- factor(trajs.sum$mark, levels = c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"))

jsub <- trajs.sum %>% filter(pos > 0e7 & pos < 999e7)
jlims <- range(jsub$exprs)
jmid <- min(jlims) + (max(jlims) - min(jlims)) / 2



pdf(file = paste0("~/data/scchic/pdfs/variance_over_pseudotime_heatmaps_zoomed_in.", Sys.Date(), ".pdf"), useDingbats = FALSE)
for (w in sort(unique(jsub$pos.round))){
  print(paste("w:", w))
  jsubsub <- jsub %>% filter(pos.round == w)
  m.lines <- ggplot(jsubsub, aes(x = pos / 1e6, y = exprs)) + geom_line() + facet_grid(lambda.bin ~ mark) + 
    theme_bw() + 
    xlab("Position (MB)") + ylab("Log2 Exprs") + 
    theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  m1 <- ggplot(jsubsub, aes(x = pos / 1e6, y = reorder(lambda.bin, desc(lambda.bin)), fill = exprs)) + 
    # geom_tile(width = 0.1, height = 1) + 
    geom_tile() + 
    theme_bw() + theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_fill_gradient(high = "darkblue", low = "gray85") + 
    scale_y_discrete(breaks = seq(0, 1, length.out = 2)) + 
    facet_wrap(~mark, ncol = 1) + 
    ggtitle(paste(jtraj, jstr, w)) + 
    ylab("Trajectory") + 
    xlab("Position (MB)")
  m1.rev <- ggplot(jsubsub, aes(x = pos / 1e6, y = reorder(lambda.bin, lambda.bin), fill = exprs)) + 
    # geom_tile(width = 0.1, height = 1) + 
    geom_tile() + 
    theme_bw() + theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_fill_gradient(high = "darkblue", low = "gray85") + 
    scale_y_discrete(breaks = seq(0, 1, length.out = 2)) + 
    facet_wrap(~mark, ncol = 1) + 
    ggtitle(paste(jtraj, jstr, w)) + 
    ylab("Trajectory") + 
    xlab("Position (MB)")
  print(m.lines)
  print(m1)
  print(m1.rev)
}

m1.all <- ggplot(jsub, aes(x = pos / 1e6, y = reorder(lambda.bin, desc(lambda.bin)), fill = exprs)) + 
  # geom_tile(width = 0.1, height = 1) + 
  geom_tile() + 
  theme_bw() + theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_gradient(high = "darkblue", low = "gray85") + 
  facet_wrap(~mark, ncol = 1) + 
  ggtitle(paste(jtraj, jstr, w)) + 
  ylab("Trajectory") + 
  xlab("Position (MB)")
m1.rev.all <- ggplot(jsub, aes(x = pos / 1e6, y = reorder(lambda.bin, lambda.bin), fill = exprs)) + 
  # geom_tile(width = 0.1, height = 1) + 
  geom_tile() + 
  theme_bw() + theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_gradient(high = "darkblue", low = "gray85") + 
  facet_wrap(~mark, ncol = 1) + 
  ggtitle(paste(jtraj, jstr, w)) + 
  ylab("Trajectory") + 
  xlab("Position (MB)")

# show mutual exclusiivity between the two repressive marks 
m.slopes <- ggplot(jsub.fits, aes(x = pos / 10^6, y = reorder(mark, desc(mark)), fill = slope)) + geom_tile() + 
  theme_bw() + theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  scale_fill_gradient2(low = "darkred", mid = "gray85", high = "darkblue", midpoint = 0, name = "Slope [A.U.]") + 
  ylab("") + 
  xlab("Position (MB)")
print(m.slopes)
# plot mutual exclusivity
m.cor <- ggpairs(jsub.wide %>% select(-pos),
                 lower = list(continuous = wrap("points", alpha = 0.2, size = 0.5))) + theme_classic() + ggtitle(jstr)
print(m.cor)

dev.off()


# w <- 40000000
# jsub.fits <- jsub %>% filter(pos.round == w) %>%


# ggplot(jsub.wide, aes(x = H3K27me3, y = H3K9me3)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# jstrs <- paste("chr", seq(21), ":", sep = "")
# 
# for (jstr in jstrs){
#   print(jstr)
# }

# Do SD over time ---------------------------------------------------------
# 
# # all trajs
# 
# # get genome wide bins: by marks 
# jmark <- "H3K4me1"
# mat.sub.gw <- GetMatSub(tm.result.lst, jmark, "", jpseudo, jfac) %>% mutate(mark = jmark)
# 
