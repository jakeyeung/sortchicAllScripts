# Jake Yeung
# Date of Creation: 2019-05-14
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/8-sharpening_over_pseudotime.R
# Do sharpening over pseudotime 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

LoadGetImputed <- function(inf){
  load(inf)  # out.objs, dat.umap.long, custom.settings. louv.settings
  imput.mat <- t(out.objs$tm.result$topics %*% out.objs$tm.result$terms)
  return(imput.mat)
}

LoadGetTmResult <- function(inf){
  load(inf)  # out.objs, dat.umap.long, custom.settings. louv.settings
  return(out.objs$tm.result)
}

source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/PlotFunctions.R")

# Load LDA objects --------------------------------------------------------

indir <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmark <- jmarks[["H3K4me1"]]
infs <- lapply(jmarks, function(jmark) list.files(indir, pattern = paste0(jmark, ".RData"), full.names = TRUE))

tm.result.lst <- lapply(infs[jmark], LoadGetTmResult)


# Load trajectories -------------------------------------------------------

inf.traj <- paste0("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_", jmark, ".Rdata")
load(inf.traj, v=T)

trajs.mark <- trajs[[jmark]]


# Make long for chromosome 15 ---------------------------------------------

jpseudo <- 0
jfac <- 10^6
jstr <- "chr15:"

mat.sub.merge <- lapply(jmark, function(jmark) GetMatSub(tm.result.lst, jmark, jstr, jpseudo, jfac) %>% mutate(mark = jmark)) %>%
  bind_rows()


# Do one trajectory -------------------------------------------------------


jtraj <- "granu"

trajs.long <- lapply(trajs, function(x) x[[jtraj]]) %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(mark = strsplit(cell, "-")[[1]][[4]]) %>%
  left_join(., mat.sub.merge) %>%
  rowwise() %>%
  mutate(lambda.bin = floor(lambda * 10) / 10)
trajs.sum <- trajs.long %>%
  group_by(lambda.bin, mark, coord, pos) %>%
  summarise(exprs = mean(exprs))

chunksize <- 2e7  # 20 MB
trajs.sum$pos.round <- as.integer(floor(trajs.sum$pos / chunksize) * chunksize)

jsub <- trajs.sum %>% filter(pos > 0e7 & pos < 999e7)
jlims <- range(jsub$exprs)
jmid <- min(jlims) + (max(jlims) - min(jlims)) / 2

m.lines <- ggplot(trajs.sum, aes(x = pos, y = exprs)) + geom_line() + facet_grid(lambda.bin ~ mark) +
  theme_bw() +
  theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

rounds <- unique(jsub$pos.round)
w <- rounds[[4]]
jsubsub <- jsub %>% filter(pos.round == w)

m1 <- ggplot(jsubsub, aes(x = pos / 1e6, y = reorder(lambda.bin, dplyr::desc(lambda.bin)), fill = exprs)) +
  geom_tile() +
  theme_bw() + theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
  scale_fill_gradientn(colours = colorRamps::matlab.like(5)) +
  scale_y_discrete(breaks = seq(0, 1, length.out = 2)) +
  facet_wrap(~mark, ncol = 1) +
  ggtitle(paste(jtraj, jstr, w)) +
  ylab("Trajectory") +
  xlab("Position (MB)")
print(m.lines)
print(m1)
