# Jake Yeung
# Date of Creation: 2019-05-26
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/make_figs/make_fig4.R
# Make figures for figure 3



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


indir <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
inf.traj <- paste0("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata")

assertthat::assert_that(dir.exists(indir))
assertthat::assert_that(file.exists(inf.traj))

pdfout <- paste0("/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory/variance_over_pseudotime_fewer_trajs.", Sys.Date(), ".pdf")

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

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmark <- jmarks[["H3K4me1"]]
infs <- lapply(jmarks, function(jmark) list.files(indir, pattern = paste0(jmark, ".RData"), full.names = TRUE))
tm.result.lst <- lapply(infs, LoadGetTmResult)


# Load trajectories -------------------------------------------------------

load(inf.traj, v=T)


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

pdf(pdfout, useDingbats = FALSE)
print(m.facet)
print(m.nofacet)
dev.off()