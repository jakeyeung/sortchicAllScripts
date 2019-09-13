# Jake Yeung
# Date of Creation: 2019-06-07
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/filter_stringent_H3K4me3/make_scChIC_signal_across_chromosome.R
# Show scChIC signal across chromosome 



library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(gganimate)

source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/AuxB6.R")
source("scripts/Rfunctions/PlotFunctions.R")


# Constants ---------------------------------------------------------------

# outdir <- "/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory"
outdir <- "/Users/yeung/data/scchic/pdfs/B6_figures/stringent_pdfs/trajectories_stringent"

# Load LDA objects --------------------------------------------------------

indir <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmark <- jmarks[["H3K4me1"]]
infs <- lapply(jmarks, function(jmark) list.files(indir, pattern = paste0(jmark, ".RData"), full.names = TRUE))
infs.stringent <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.H3K4me3.stringent_filter.RData"

infs$H3K4me3 <- infs.stringent

tm.result.lst <- lapply(infs, LoadGetTmResult)


# Load trajectories -------------------------------------------------------

inf.traj.stringent <- "/Users/yeung/data/scchic/robjs/B6_objs_stringent/traj_objs_H3K4me3.stringent.Rdata"
load(inf.traj.stringent, v=T)
dat.umap.long.trajs.stringent <- dat.umap.long.trajs$H3K4me3
trajs.stringent <- trajs$H3K4me3
trajs.objs.stringent <- trajs.objs$H3K4me3

inf.traj <- "/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata"
load(inf.traj, v=T)

dat.umap.long.trajs$H3K4me3 <- dat.umap.long.trajs.stringent
trajs$H3K4me3 <- trajs.stringent
trajs.objs$H3K4me3 <- trajs.objs.stringent


# Make long for chromosome 15 ---------------------------------------------

jpseudo <- 0
jfac <- 10^6
jstr <- "chr15:"

mat.sub.merge <- lapply(jmarks, function(jmark) GetMatSub(tm.result.lst, jmark, jstr, jpseudo, jfac) %>% mutate(mark = jmark)) %>%
  bind_rows()


# Do one trajectory -------------------------------------------------------

chunksize <- 2e7  # 20 MB
trajnames <- c("lymphoid", "granu", "eryth", "mega")
# jtraj <- "eryth"

# do one trajectory, two marks 
# jtraj <- "granu"
# jtraj <- "lymphoid"
jtraj <- "eryth"



print(jtraj)
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

trajs.sum$mark <- factor(as.character(trajs.sum$mark), levels = jmarks)

trajs.sum$pos.round <- as.integer(floor(trajs.sum$pos / chunksize) * chunksize)

jlims <- range(trajs.sum$exprs)
jmid <- min(jlims) + (max(jlims) - min(jlims)) / 2

rounds <- unique(trajs.sum$pos.round)
w <- rounds[[4]]


pdf(paste0("/Users/yeung/data/scchic/pdfs/B6_figures/stringent_pdfs/variance_over_trajectory_stringent/trajectory_lines_over_pseudotime.", jtraj, ".", Sys.Date(), ".pdf"), useDingbats = FALSE)
m.lines <- ggplot(subset(trajs.sum, pos.round == w), aes(x = pos / 10^6, y = exprs, color = exprs)) + geom_line() + 
  facet_grid(lambda.bin ~ mark) +
  theme_bw(20) +
  theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_blank(), legend.position = "bottom") + 
  ggtitle(paste(jtraj, jstr)) + 
  scale_color_gradientn(colours = colorRamps::matlab.like(5)) + 
  xlab("Position [MB]") + ylab(expression(log[2]~"(signal)"))
print(m.lines)

# do for separate marks
for (jmark in jmarks){
  m.lines.mark <- ggplot(subset(trajs.sum, pos.round == w & mark == jmark), aes(x = pos / 10^6, y = exprs, color = exprs)) + geom_line() + 
    facet_grid(lambda.bin ~ mark) +
    theme_bw(20) +
    theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_blank(), legend.position = "bottom") + 
    ggtitle(paste(jtraj, jstr)) + 
    scale_color_gradientn(colours = colorRamps::matlab.like(5)) + 
    xlab("Position [MB]") + ylab(expression(log[2]~"(signal)"))
  print(m.lines.mark)
}

# do for separate marks, lambda = 0, lambda = 1
for (jmark in jmarks){
  m.lines.mark.pseudo <- ggplot(subset(trajs.sum, pos.round == w & lambda.bin %in% c(0, 1) & mark == jmark), aes(x = pos / 10^6, y = exprs, color = exprs)) + geom_line() + facet_grid(lambda.bin ~ mark) +
    theme_bw(20) +
    theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_blank(), legend.position = "bottom") + 
    # ggtitle(paste(jtraj, jstr)) + 
    scale_color_gradientn(colours = colorRamps::matlab.like(5)) + 
    xlab("Position [MB]") + ylab(expression(log[2]~"(signal)"))
  print(m.lines.mark.pseudo)
}
dev.off()

# marks.sub.lst <- list(c("H3K4me1", "H3K9me3"), c("H3K27me3", "H3K9me3"), c("H3K4me1", "H3K9me3"), c("H3K4me1", "H3K27me3"), c("H3K4me1", "H3K27me3"), c("H3K4me1", "H3K4me3"))
marks.sub.lst <- list(c("H3K4me1", "H3K9me3"), c("H3K27me3", "H3K9me3"), c("H3K4me1", "H3K4me3"))
for (marks.sub in marks.sub.lst){
  print(marks.sub)
  marks.str <- paste(marks.sub, collapse = "_")
  jsub <- subset(trajs.sum, mark %in% marks.sub)
  # print pseudo 0 and pseudo 1 as anchor points, then do the gif
  
  # gifpath <- paste0("/Users/yeung/data/scchic/pdfs/B6_figures/for_presentation/sharpening.", Sys.Date(), ".", marks.str, ".gif")
  gifpath <- paste0("/Users/yeung/data/scchic/pdfs/B6_figures/stringent_pdfs/variance_over_trajectory_stringent/sharpening.", jtraj, ".", Sys.Date(), ".", marks.str, ".stringent.gif")
  assertthat::assert_that(dir.exists(dirname(gifpath)))
  # make animaition 
  anim <- ggplot(subset(jsub, pos.round == w), aes(x = pos / 10^6, y = exprs, color = exprs)) + geom_line() + facet_wrap(~mark, ncol = 1) + 
    theme_bw() +
    theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    labs(title = 'Pseudotime: {closest_state}') +    
    scale_color_gradientn(colours = colorRamps::matlab.like(5), name = "Log Signal") + 
    xlab("Position [MB]") + ylab(expression(log[2]~"(signal)")) + 
    transition_states(
      lambda.bin, transition_length = 1, state_length = 0
    ) + 
    ease_aes('linear')
  animate(anim, fps = 12, start_pause = 25, end_pause = 0, rewind = FALSE, nframes = 100, renderer = gifski_renderer(gifpath))
}
