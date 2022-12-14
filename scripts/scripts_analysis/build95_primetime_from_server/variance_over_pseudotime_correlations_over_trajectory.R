# Jake Yeung
# Date of Creation: 2019-05-03
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/variance_over_pseudotime_heatmaps_zoomed_in.R
# Zoom in on genomic regions 


rm(list=ls())

tstart <- Sys.time() 

library(data.table)
library(dplyr)
library(ggplot2)
library(JFuncs)
library(tidyr)
library(GGally)
library(gridExtra)

library(grDevices)

library(GGally)

library(parallel)

source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/TrajFunctions.R")





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


GetTrajSum <- function(tm.result.lst, trajs.mixed, jmarks, jstr, jpseudo, jfac, jtraj){
  mat.sub.merge <- lapply(jmarks, function(jmark) GetMatSub(tm.result.lst, jmark, jstr, jpseudo, jfac) %>% mutate(mark = jmark)) %>% 
    bind_rows()
  trajs.long <- lapply(trajs.mixed, function(x) x[[jtraj]]) %>%
    bind_rows() %>%
    rowwise() %>%
    mutate(mark = strsplit(cell, "_")[[1]][[2]]) %>%
    left_join(., mat.sub.merge) %>%
    rowwise() %>%
    mutate(lambda.bin = floor(lambda * 10) / 10)
  trajs.sum <- trajs.long %>%
    group_by(lambda.bin, mark, coord, pos) %>%
    summarise(exprs = mean(exprs)) %>%
    mutate(chromo = jstr)
  return(trajs.sum)
}

FitSlope <- function(dat.sub){
  # fit exprs to trajectory to find slope
  fit <- lm(formula = exprs ~ lambda.bin, data = dat.sub)
  slope <- fit$coefficients[['lambda.bin']]
  int <- fit$coefficients[['(Intercept)']]
  pval <- summary(fit)$coefficients["lambda.bin", "Pr(>|t|)"]
  return(data.frame(slope = slope, int = int, pval = pval))
}

# Load data  --------------------------------------------------------------

inf.trajs <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/primeitme_objects/trajectory_from_spring_2019-04-15.RData"
assertthat::assert_that(file.exists(inf.trajs))
load(inf.trajs, v=T)

inf.dat <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/primeitme_objects/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"
assertthat::assert_that(file.exists(inf.dat))
load(inf.dat, v=T)

trajs.mixed.lst <- GetTrajMixed(inf.trajs, inf.dat)
trajs.spring <- trajs.mixed.lst$trajs.mixed
trajs.mixed <- trajs.mixed.lst$trajs.mixed

# Get constants -----------------------------------------------------------


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


jstr <- ""

jthres <- 0.05

# Explore the data  -------------------------------------------------------

jmark <- "H3K4me1"
jmark <- "H3K27me3"
jmark <- "H3K9me3"

jtraj <- "granu"

# Prepare dat -------------------------------------------------------------

jpseudo <- 0
jfac <- 10^6
jtraj <- "granu"
chunksize <- 2e7  # 20 MB

jstrs <- paste("chr", seq(21), ":", sep = "")
# jstrs <- c("chr15:")

# this takes long time! Can we somehow speed it up?
trajs.sum <- mclapply(jstrs[1:11], function(jstr) GetTrajSum(tm.result.lst = tm.result.lst, trajs.mixed = trajs.mixed, jmarks, jstr = jstr, jfac = jfac, jpseudo = jpseudo, jtraj = jtraj), mc.cores = 11) %>%
  bind_rows()

outf <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/robjs/trajs_sum_variance_from_server.rds"
if (!file.exists(outf)){
  saveRDS(trajs.sum, file = outf)
}

# trajs.sum <- trajs.sum %>% mutate(chromo = jstrs[[1]])

# Do heatmap of granulocytes for the 4 marks ------------------------------

jsub.fits <- trajs.sum %>%
  group_by(pos, mark, chromo) %>%
  do(FitSlope(.))

jsub.fits$mark <- factor(jsub.fits$mark, levels = jmarks)

jsub.wide <- subset(jsub.fits, select = c(pos, mark, slope, chromo)) %>% spread(key = mark, value = slope) %>% ungroup()

trajs.sum$pos.round <- as.integer(floor(trajs.sum$pos / 2e7) * 2e7)
trajs.sum$pos.round2 <- as.integer(floor(trajs.sum$pos / 5e6) * 5e6)
trajs.sum$mark <- factor(trajs.sum$mark, levels = c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"))

jsub <- trajs.sum %>% filter(pos > 0e7 & pos < 999e7)
jlims <- range(jsub$exprs)
jmid <- min(jlims) + (max(jlims) - min(jlims)) / 2

pdf(file = paste0("/hpc/hub_oudenaarden/jyeung/data/scChiC/pdfs/variance_over_pseudotime_slopes_and_correlations_genomewide.", Sys.Date(), ".pdf"), useDingbats = FALSE)
  # show mutual exclusiivity between the two repressive marks 
  for (jchromo in jstrs){
    m.slopes <- ggplot(jsub.fits %>% filter(chromo == jchromo), aes(x = pos / 10^6, y = reorder(mark, desc(mark)), fill = slope)) + geom_tile() + 
      theme_bw() + theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
      scale_fill_gradient2(low = "darkred", mid = "gray85", high = "darkblue", midpoint = 0, name = "Slope [A.U.]") + 
      ylab("") + 
      xlab("Position (MB)") + 
      ggtitle(jchromo)
      # facet_wrap(~chromo)
    print(m.slopes)
  }
  # plot mutual exclusivity
  m.cor <- ggpairs(jsub.wide %>% select(-pos, -chromo),
                   lower = list(continuous = HexPlot)) + theme_classic() + ggtitle("GenomeWide")
                   # lower = list(continuous = wrap("points", alpha = 0.2, size = 0.5))) + theme_classic() + ggtitle(jstr)
  print(m.cor)
dev.off()


print(Sys.time() - tstart)
