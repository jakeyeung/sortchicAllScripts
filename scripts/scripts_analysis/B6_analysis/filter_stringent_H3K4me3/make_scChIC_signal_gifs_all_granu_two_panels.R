# Jake Yeung
# Date of Creation: 2019-09-22
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/filter_stringent_H3K4me3/make_scChIC_signal_gifs_all_granu_two_panels.R
# description


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(gganimate)

library(magick)

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
jtraj <- "granu"
# jtraj <- "lymphoid"
# jtraj <- "eryth"


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

trajs.long.nopos <- lapply(trajs, function(x) x[[jtraj]]) %>%
  bind_rows() %>%
  rowwise() %>%
  mutate(mark = strsplit(cell, "-")[[1]][[4]]) %>%
  rowwise() %>%
  mutate(lambda.bin = floor(lambda * 10) / 10)

pos.sum <- trajs.long %>%
  group_by(lambda.bin, mark) %>%
  summarise(umap1 = mean(umap1),
            umap2 = mean(umap2))

ggplot(pos.sum, aes(x = umap1, y = umap2, color = lambda.bin)) + geom_point() + scale_color_viridis_c() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

trajs.sum$mark <- factor(as.character(trajs.sum$mark), levels = jmarks)

trajs.sum$pos.round <- as.integer(floor(trajs.sum$pos / chunksize) * chunksize)

jlims <- range(trajs.sum$exprs)
jmid <- min(jlims) + (max(jlims) - min(jlims)) / 2

rounds <- unique(trajs.sum$pos.round)
w <- rounds[[4]]

jlambdas <- sort(unique(trajs.sum$lambda.bin))
names(jlambdas) <- jlambdas
dat.umap.long.trajs.lambda <- lapply(jlambdas, function(jlamb){
  jtmp <- bind_rows(dat.umap.long.trajs)
  jtmp$lambda <- jlamb
  return(subset(jtmp, select = c(cell, umap1, umap2, mark, lambda)))
}) %>%
  bind_rows()
dat.umap.long.trajs.lambda$mark <- factor(as.character(dat.umap.long.trajs.lambda$mark), levels = jmarks)

trajs.long.nopos.lambda <- lapply(jlambdas, function(jlamb){
  jtmp <- trajs.long.nopos
  jtmp$lambda <- jlamb
  return(subset(jtmp, select = c(cell, umap1, umap2, mark, lambda)))
}) %>%
  bind_rows()
trajs.long.nopos.lambda$mark <- factor(as.character(trajs.long.nopos.lambda$mark), levels = jmarks)

pdf(paste0("/Users/yeung/data/scchic/pdfs/B6_figures/stringent_pdfs/variance_over_trajectory_stringent/trajectory_lines_over_pseudotime.", jtraj, ".", Sys.Date(), ".pdf"), useDingbats = FALSE)
m.lines <- ggplot(subset(trajs.sum, pos.round == w), aes(x = pos / 10^6, y = exprs, color = exprs)) + geom_line() + 
  facet_grid(lambda.bin ~ mark) +
  theme_bw(20) +
  theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_blank(), legend.position = "bottom") + 
  ggtitle(paste(jtraj, jstr)) + 
  scale_color_gradientn(colours = colorRamps::matlab.like(5)) + 
  xlab("Position [MB]") + ylab(expression(log[2]~"(signal)"))
print(m.lines)



# dat.umap.long.trajs.lambda <- 
m.dots <- ggplot(dat.umap.long.trajs.lambda, aes(x = umap1, y = umap2)) + 
  geom_point(alpha = 0.1, color = "gray80") + 
  # geom_path(data = trajs.long.nopos.lambda, inherit.aes = FALSE, aes(x = umap1, y = umap2, color = "green"), size = 1) + 
  geom_point(mapping = aes(x = umap1, y = umap2), data = pos.sum, size = 2, color = 'blue') + 
  facet_grid(lambda.bin ~ mark) +
  theme_bw(20) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_blank(), legend.position = "bottom") 
print(m.dots)

  # scale_color_viridis_c()
  # scale_color_gradientn(colours = colorRamps::matlab.like(5)) 


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
  # jsub.pos <- subset(pos.sum, mark %in% marks.sub)
  # print pseudo 0 and pseudo 1 as anchor points, then do the gif
  
  # gifpath <- paste0("/Users/yeung/data/scchic/pdfs/B6_figures/for_presentation/sharpening.", Sys.Date(), ".", marks.str, ".gif")
  gifpath <- paste0("/Users/yeung/data/scchic/pdfs/B6_figures/stringent_pdfs/variance_over_trajectory_stringent/sharpening.", jtraj, ".", Sys.Date(), ".", marks.str, ".stringent.gif")
  assertthat::assert_that(dir.exists(dirname(gifpath)))
  # make animaition 
  anim <- ggplot(subset(jsub, pos.round == w), aes(x = pos / 10^6, y = exprs, color = exprs)) + geom_line() + facet_wrap(~mark, ncol = 1) + 
    theme_bw() +
    theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    labs(title = 'Pseudotime: {closest_state}') +    
    scale_color_gradientn(colours = colorRamps::matlab.like(5), name = "Log Signal") + 
    xlab("Position [MB]") + ylab(expression(log[2]~"(signal)")) + 
    transition_states(
      lambda.bin, transition_length = 1, state_length = 0
    ) + 
    ease_aes('linear')
  # animate(anim, fps = 12, start_pause = 25, end_pause = 0, rewind = FALSE, nframes = 100, renderer = gifski_renderer(gifpath))
  anim.out1 <- animate(anim, fps = 12, start_pause = 25, end_pause = 0, rewind = FALSE, nframes = 100, width = 240 * 4, height = 240 * 2)
  
  
  # m.dots <- ggplot(dat.umap.long.trajs.lambda, aes(x = umap1, y = umap2)) + 
  #   geom_point(alpha = 0.1, color = "gray80") + 
  #   geom_point(mapping = aes(x = umap1, y = umap2), data = pos.sum, size = 2, color = 'blue') + 
  #   facet_grid(lambda.bin ~ mark) +
  #   theme_bw(20) +
  #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_blank(), legend.position = "bottom") 
  # print(m.dots)
  
  anim2 <- ggplot(subset(dat.umap.long.trajs.lambda, mark %in% marks.sub), aes(x = umap1, y = umap2)) + facet_wrap(~mark, ncol = 1, scales = "free") + 
    geom_point(alpha = 0.1, color = "gray80") + 
    geom_point(mapping = aes(x = umap1, y = umap2), data = subset(pos.sum, mark %in% marks.sub), size = 2, color = 'blue') + 
    theme_minimal() +
    xlab("") + ylab("") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_blank(), legend.position = "bottom", axis.ticks = element_blank(), axis.text = element_blank())  + 
    transition_states(
      lambda.bin, transition_length = 1, state_length = 0
    ) + 
    ease_aes('linear')
  
  anim.out2 <- animate(anim2, fps = 12, start_pause = 25, end_pause = 0, rewind = FALSE, nframes = 100, width = 240, height = 240 * 2)
  
  d<-image_blank(240*3, 240 * 2)
  the_frame<-d
  for(i in 2:100){
    the_frame<-c(the_frame,d)
  }
  
  a_mgif<-image_read(anim.out1)
  b_mgif<-image_read(anim.out2)
  
  new_gif<-image_append(c(a_mgif[1], b_mgif[1]))
  for(i in 2:100){
    combined <- image_append(c(a_mgif[i], b_mgif[i]))
    new_gif<-c(new_gif,combined)
  }
  magick::image_write(new_gif, path=gifpath)
  # do the other animation
  
}

