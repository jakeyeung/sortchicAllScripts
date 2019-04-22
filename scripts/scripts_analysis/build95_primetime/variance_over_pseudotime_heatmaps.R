# Jake Yeung
# Date of Creation: 2019-04-21
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/variance_over_pseudotime_heatmaps.R
# Plot variance along pseudotime with heatmaps 


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

source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/Aux.R")



# Load data  --------------------------------------------------------------

# inf.trajs <- "/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-09.RData"
inf.trajs <- "/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-11.RData"
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
jfac <- 10^7


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

ggplot(trajs.sum, aes(x = pos, y = exprs)) + geom_line() + facet_grid(lambda.bin ~ mark) + 
  theme_bw() + 
  theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jsub <- trajs.sum %>% filter(pos > 0e7 & pos < 999e7)
jlims <- range(jsub$exprs)
jmid <- min(jlims) + (max(jlims) - min(jlims)) / 2
m1 <- ggplot(jsub, aes(x = pos / 1e6, y = lambda.bin, fill = exprs)) + 
  # geom_tile(width = 0.1, height = 1) + 
  geom_tile() + 
  theme_bw() + theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_fill_gradient2(low = "gray85", mid = "gray50", high = scales::muted("darkblue"), midpoint = jmid, lim = jlims, name = "log2(exprs)") + 
  facet_wrap(~mark, ncol = 1) + 
  ggtitle(paste(jtraj, jstr)) + 
  ylab("Trajectory") + 
  xlab("Position (MB)")
print(m1)

 
# # do heatmap
# 
# jmat <- trajs.sum %>% 
#   filter(mark == jmark) %>% 
#   ungroup() %>%
#   dplyr::select(pos, lambda.bin, exprs) %>% 
#   spread(key = pos, value = exprs) 
# 
# 
# 
# 
