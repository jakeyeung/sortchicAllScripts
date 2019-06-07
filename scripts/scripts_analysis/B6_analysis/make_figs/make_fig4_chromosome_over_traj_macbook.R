# Jake Yeung
# Date of Creation: 2019-06-03
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/make_figs/make_fig4_chromosome_over_traj_macbook.R
# Macbook 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(GGally)
library(colorRamps)
# library(ggrastr)
library(parallel)

source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/PlotFunctions.R")


# indir <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/bamlist_for_merging_build95_B6"
# inf.traj <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/from_macbook/B6_objs/traj_objs_all_marks.Rdata"
indir <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
inf.traj <- paste0("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata")

assertthat::assert_that(dir.exists(indir))
assertthat::assert_that(file.exists(inf.traj))

# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) == 1){
#   outdir <- args[[1]]
#   assertthat::assert_that(dir.exists(outdir))
# } else {
#   outdir <- "/tmp"
# }
outdir <- "/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory"


# Constants ---------------------------------------------------------------

# outdir <- "/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory"

# Load LDA objects --------------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmark <- jmarks[["H3K4me1"]]
infs <- lapply(jmarks, function(jmark) list.files(indir, pattern = paste0(jmark, ".RData"), full.names = TRUE))

tm.result.lst <- lapply(infs, LoadGetTmResult)


# Load trajectories -------------------------------------------------------

load(inf.traj, v=T)


# Make long for chromosome 15 ---------------------------------------------

jpseudo <- 0
jfac <- 10^6
jstr <- "chr15:"
nbins <- 50
pdfout <- file.path(outdir, paste0("Fig4_SharpeningOverTraj.", nbins, ".pdf"))

mat.sub.merge <- lapply(jmarks, function(jmark) GetMatSub(tm.result.lst, jmark, jstr, jpseudo, jfac) %>% mutate(mark = jmark)) %>%
  bind_rows()


# Do one trajectory -------------------------------------------------------

chunksize <- 2e7  # 20 MB
trajnames <- c("lymphoid", "granu", "eryth", "mega")
# jtraj <- "eryth"


# for (jtraj in trajnames){
mclapply(trajnames, function(jtraj){
  print(jtraj)
  outpdf <- file.path(outdir, paste0("heatmap_exprs_over_traj_", jtraj, ".", Sys.Date(), ".pdf"))
  if (file.exists(outpdf)){
    print(paste("Skipping trajectory:", jtraj))
    next
  }
  pdf(outpdf, useDingbats = FALSE)
  trajs.long <- lapply(trajs, function(x) x[[jtraj]]) %>%
    bind_rows() %>%
    rowwise() %>%
    mutate(mark = strsplit(cell, "-")[[1]][[4]]) %>%
    left_join(., mat.sub.merge) %>%
    rowwise() %>%
    mutate(lambda.bin = floor(lambda * nbins) / nbins)
  trajs.sum <- trajs.long %>%
    group_by(lambda.bin, mark, coord, pos) %>%
    summarise(exprs = mean(exprs))
  
  trajs.sum$mark <- factor(as.character(trajs.sum$mark), levels = jmarks)
  
  trajs.sum$pos.round <- as.integer(floor(trajs.sum$pos / chunksize) * chunksize)
  
  jsub <- trajs.sum
  jlims <- range(jsub$exprs)
  jmid <- min(jlims) + (max(jlims) - min(jlims)) / 2
  
  m.lines <- ggplot(trajs.sum, aes(x = pos, y = exprs)) + geom_line() + facet_grid(lambda.bin ~ mark) +
    theme_bw() +
    theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(paste(jtraj, jstr))
  print(m.lines)
  
  rounds <- unique(jsub$pos.round)
  # w <- rounds[[4]]
  for (w in rounds){
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
    print(m1)
    m1.rev <- ggplot(jsubsub, aes(x = pos / 1e6, y = reorder(lambda.bin, lambda.bin), fill = exprs)) +
      geom_tile() +
      theme_bw() + theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
      scale_fill_gradientn(colours = colorRamps::matlab.like(5)) +
      scale_y_discrete(breaks = seq(0, 1, length.out = 2)) +
      facet_wrap(~mark, ncol = 1) +
      ggtitle(paste(jtraj, jstr, w)) +
      ylab("Trajectory") +
      xlab("Position (MB)")
    print(m1.rev)
  }
  dev.off()
# }
}, mc.cores = length(trajnames))

# 
# # Do slope calculations  --------------------------------------------------
# 
# # jtraj <- "granu"
# 
# # for (jtraj in trajnames){
# mclapply(trajnames, function(jtraj){
#   print(jtraj)
#   outpdf <- file.path(outdir, paste0("heatmap_exprs_over_traj_slopes_", jtraj, ".", Sys.Date(), ".pdf"))
#   pdf(outpdf, useDingbats = FALSE)
#   trajs.long2 <- lapply(trajs, function(x) x[[jtraj]]) %>%
#     bind_rows() %>%
#     rowwise() %>%
#     mutate(mark = strsplit(cell, "-")[[1]][[4]]) %>%
#     left_join(., mat.sub.merge) %>%
#     rowwise() %>%
#     mutate(lambda.bin = floor(lambda * 10) / 10)
#   trajs.sum <- trajs.long2 %>%
#     group_by(lambda.bin, mark, coord, pos) %>%
#     summarise(exprs = mean(exprs))
#   
#   jsub.fits <- trajs.sum %>%
#     group_by(pos, mark) %>%
#     do(FitSlope(.))
#   
#   jsub.fits$mark <- factor(jsub.fits$mark, levels = c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"))
#   jsub.ref <- subset(jsub.fits, mark == "H3K9me3") %>% ungroup() %>% dplyr::rename(slope.H3K9me3 = slope) %>% dplyr::select(pos, slope.H3K9me3)
#   jsub.compare <- subset(jsub.fits, mark != "H3K9me3")
#   jsub.compare <- left_join(jsub.compare, jsub.ref)
#   
#   jsub.wide <- subset(jsub.fits, select = c(pos, mark, slope)) %>% spread(key = mark, value = slope) %>% ungroup()
#   
#   
#   # show mutual exclusiivity between the two repressive marks
#   m.slopes <- ggplot(jsub.fits, aes(x = pos / 10^6, y = reorder(mark, dplyr::desc(mark)), fill = slope)) + geom_tile() +
#     theme_bw() + theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
#     # scale_fill_gradient2(low = "darkred", mid = "gray85", high = "darkblue", midpoint = 0, name = "Slope [A.U.]") +
#     # scale_fill_gradient2(low = "darkred", mid = "gray85", high = "darkblue", midpoint = 0, name = "Slope [A.U.]") +
#     scale_fill_gradientn(colours = colorRamps::matlab.like(5), name = "Slope [A.U.]") +
#     ylab("") +
#     xlab("Position (MB)")
#   print(m.slopes)
#   # plot mutual exclusivity
#   m.cor <- ggpairs(jsub.wide %>% dplyr::select(-pos),
#                    lower = list(continuous = wrap("points", alpha = 0.2, size = 0.5))) + theme_classic() + ggtitle(jstr)
#   print(m.cor)
#   # do H3K27me3 vs H3K9me3 only
#   m.cor2 <- ggpairs(subset(jsub.wide %>% dplyr::select(H3K9me3, H3K27me3, -pos)), lower = list(continuous = wrap("points", alpha = 0.25, size = 1.5))) + theme_classic(20) + ggtitle(jstr)
#   print(m.cor2)
#   
#   m <- ggplot(jsub.compare, aes(y = slope.H3K9me3, x = slope)) +
#     # geom_point_rast(alpha = 0.2, size = 6) +
#     geom_point(alpha = 0.2) +
#     geom_density2d() + facet_wrap(~mark, scales = "free_x") +
#     theme_bw(18) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#     xlab("Slope") + ylab("Slope (H3K9me3)")
#   print(m)
#   
#   dev.off()
# # }
# }, mc.cores = length(trajnames))
