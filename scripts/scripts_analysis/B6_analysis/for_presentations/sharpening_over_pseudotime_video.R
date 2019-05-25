# Jake Yeung
# Date of Creation: 2019-05-19
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/for_presentations/sharpening_over_pseudotime_video.R
# Sharp over pseudotime 


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

outdir <- "/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory"

# Load LDA objects --------------------------------------------------------

indir <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmark <- jmarks[["H3K4me1"]]
infs <- lapply(jmarks, function(jmark) list.files(indir, pattern = paste0(jmark, ".RData"), full.names = TRUE))

tm.result.lst <- lapply(infs, LoadGetTmResult)


# Load trajectories -------------------------------------------------------

inf.traj <- paste0("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata")
load(inf.traj, v=T)




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
jtraj <- "granu"



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
  
  m.lines <- ggplot(subset(trajs.sum, pos.round == w), aes(x = pos / 10^6, y = exprs, color = exprs)) + geom_line() + facet_grid(lambda.bin ~ mark) +
    theme_bw() +
    theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(paste(jtraj, jstr)) + 
    scale_color_gradientn(colours = colorRamps::matlab.like(5)) + 
    xlab("Position (MB)") + ylab("Log2 scChIC Signal")
  print(m.lines)


marks.sub.lst <- list(c("H3K4me1", "H3K9me3"), c("H3K27me3", "H3K9me3"), c("H3K4me1", "H3K9me3"), c("H3K4me1", "H3K27me3"), c("H3K4me1", "H3K27me3")) 
for (marks.sub in marks.sub.lst){
  print(marks.sub)
  marks.str <- paste(marks.sub, collapse = "_")
  gifpath <- paste0("/Users/yeung/data/scchic/pdfs/B6_figures/for_presentation/sharpening.", Sys.Date(), ".", marks.str, ".gif")
  assertthat::assert_that(dir.exists(dirname(gifpath)))
  
  # make animaition 
  jsub <- subset(trajs.sum, mark %in% marks.sub)
  anim <- ggplot(subset(jsub, pos.round == w), aes(x = pos / 10^6, y = exprs, color = exprs)) + geom_line() + facet_wrap(~mark, ncol = 1) + 
    theme_bw() +
    theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    labs(title = 'Pseudotime: {closest_state}') +    
    scale_color_gradientn(colours = colorRamps::matlab.like(5), name = "Log Signal") + 
    xlab("Position (MB)") + ylab("Log2 scChIC Signal") + 
    transition_states(
      lambda.bin, transition_length = 1, state_length = 0
    ) + 
    ease_aes('linear')
  animate(anim, fps = 12, start_pause = 25, end_pause = 0, rewind = FALSE, nframes = 100, renderer = gifski_renderer(gifpath))
}

# marks.sub <- c("H3K4me1", "H3K9me3")
