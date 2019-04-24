# Jake Yeung
# Date of Creation: 2019-04-23
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/variance_over_pseudotime_standarderror.R
# StandardError

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
source("scripts/Rfunctions/PlotFunctions.R")


# Get constants -----------------------------------------------------------

make.plots <- TRUE

colhash <- GetTrajColors(as.hash = TRUE)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")
colvec <- c("gray85", "gray50", "blue")  
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

jpseudo <- 0
jfac <- 10^6


# Load data  --------------------------------------------------------------

# trajectories from 

inf.trajs <- "/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-15.RData"
assertthat::assert_that(file.exists(inf.trajs))
load(inf.trajs, v=T)

head(trajs.spring[[4]][[1]])

inf.dat <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"
assertthat::assert_that(file.exists(inf.dat))
load(inf.dat, v=T)

rm(count.mat.lst)

# get cell counts
inf.countmat <- "/Users/yeung/data/scchic/robjs/count_mat_lst_unnorm_cellmin_550_cellmax_500000.binarize.FALSE.2019-04-23.rds"
count.mat.lst <- readRDS(inf.countmat)
count.cells <- lapply(jmarks, function(jmark){
  x <- count.mat.lst[[jmark]]
  jcounts <- Matrix::colSums(x)
  jcounts.dat <- data.frame(cell = names(jcounts), 
                            counts = jcounts,
                            countslog = log10(jcounts), 
                            mark = jmark)
}) %>% bind_rows()
rm(count.mat.lst)

trajs.mixed <- c(trajs["H3K4me1"], trajs.spring[c("H3K4me3", "H3K27me3", "H3K9me3")])

# rename bcell trajectory as lymphoid
names(trajs.mixed$H3K4me1)[which(names(trajs.mixed$H3K4me1) == "lymphoid")] <- "lymphod.naive"
names(trajs.mixed$H3K4me1)[which(names(trajs.mixed$H3K4me1) == "bcell")] <- "lymphoid"



# Plot UMAP ---------------------------------------------------------------

head(dat.trajs.long)
head(dat.umap.long.trajs[[1]])

# plot no colors
jsize <- 2
jcol <- "gray70"

dat.umap.mixed <- list(H3K4me1 = dat.umap.long.trajs[["H3K4me1"]], 
                    H3K4me3 = dat.trajs.long %>% filter(mark == "H3K4me3"),
                    H3K27me3 = dat.trajs.long %>% filter(mark == "H3K27me3"), 
                    H3K9me3 = dat.trajs.long %>% filter(mark == "H3K9me3"))
# add cell counts
dat.umap.mixed <- lapply(dat.umap.mixed, function(x){
  left_join(x, count.cells)
})

jrange <- c(2.8, 6)
if (make.plots){
  pdf(paste0("~/data/scchic/pdfs/variance_over_pseudotime_plots_primetime.Figure4.umaps.", Sys.Date(), ".pdf"), useDingbats = FALSE)
}


jleg.name <- "log10(# Cuts)"
m.h3k4me1 <- PlotXYWithColor(dat.umap.mixed[["H3K4me1"]], xvar = "umap1", yvar = "umap2", cname = "countslog", jsize = jsize, leg.name = jleg.name, jjrange = jrange)
m.h3k4me3 <- PlotXYWithColor(dat.umap.mixed[["H3K4me3"]], xvar = "X1", yvar = "X2", cname = "countslog", jsize = jsize, leg.name = jleg.name, jjrange = jrange)
m.h3k27me3 <- PlotXYWithColor(dat.umap.mixed[["H3K27me3"]], xvar = "X1", yvar = "X2", cname = "countslog", jsize = jsize, leg.name = jleg.name, jjrange = jrange)
m.h3k9me3 <- PlotXYWithColor(dat.umap.mixed[["H3K9me3"]], xvar = "X1", yvar = "X2", cname = "countslog", jsize = jsize, leg.name = jleg.name, jjrange = jrange)

m.h3k4me1.auto <- PlotXYWithColor(dat.umap.mixed[["H3K4me1"]], xvar = "umap1", yvar = "umap2", cname = "countslog", jsize = jsize, leg.name = jleg.name, jjrange = "auto")
m.h3k4me3.auto <- PlotXYWithColor(dat.umap.mixed[["H3K4me3"]], xvar = "X1", yvar = "X2", cname = "countslog", jsize = jsize, leg.name = jleg.name, jjrange = "auto")
m.h3k27me3.auto <- PlotXYWithColor(dat.umap.mixed[["H3K27me3"]], xvar = "X1", yvar = "X2", cname = "countslog", jsize = jsize, leg.name = jleg.name, jjrange = "auto")
m.h3k9me3.auto <- PlotXYWithColor(dat.umap.mixed[["H3K9me3"]], xvar = "X1", yvar = "X2", cname = "countslog", jsize = jsize, leg.name = jleg.name, jjrange = "auto")

m2.h3k4me1 <- PlotXYNoColor(dat.umap.mixed[["H3K4me1"]], xvar = "umap1", yvar = "umap2", jcol = jcol, jsize = jsize)
m2.h3k4me3 <- PlotXYNoColor(dat.umap.mixed[["H3K4me3"]], xvar = "X1", yvar = "X2", jcol = jcol, jsize = jsize)
m2.h3k27me3 <- PlotXYNoColor(dat.umap.mixed[["H3K27me3"]], xvar = "X1", yvar = "X2", jcol = jcol, jsize = jsize)
m2.h3k9me3 <- PlotXYNoColor(dat.umap.mixed[["H3K9me3"]], xvar = "X1", yvar = "X2", jcol = jcol, jsize = jsize)

# plot multi
multiplot(m.h3k4me1, m.h3k4me3, m.h3k27me3, m.h3k9me3, cols = 4)
multiplot(m.h3k4me1.auto, m.h3k4me3.auto, m.h3k27me3.auto, m.h3k9me3.auto, cols = 4)
multiplot(m2.h3k4me1, m2.h3k4me3, m2.h3k27me3, m2.h3k9me3, cols = 4)

if (make.plots){
  dev.off()
}

# Plot heatmkap -----------------------------------------------------------

if (make.plots){
  pdf(paste0("~/data/scchic/pdfs/variance_over_pseudotime_plots_primetime.Figure4.Trajs.", Sys.Date(), ".pdf"), useDingbats = FALSE)
}

mat.sub.merge <- lapply(jmarks, function(jmark) GetMatSub(tm.result.lst, jmark, jstr, jpseudo, jfac) %>% mutate(mark = jmark)) %>% 
  bind_rows()

# jtrajs <- c("granu")
jtrajs <- c("granu", "lymphoid", "eryth")

for (jtraj in jtrajs){
  print(jtraj)
  trajs.sum <- lapply(trajs.mixed, function(x) x[[jtraj]]) %>%
    bind_rows() %>%
    rowwise() %>%
    mutate(mark = strsplit(cell, "_")[[1]][[2]]) %>%
    left_join(., mat.sub.merge) %>%
    rowwise() %>%
    mutate(lambda.bin = floor(lambda * 10) / 10) %>%
    group_by(lambda.bin, mark, coord, pos) %>%
    summarise(exprs = mean(exprs))
  
  jsub <- trajs.sum %>% filter(pos > 0e7 & pos < 999e7)
  jsub$mark <- factor(jsub$mark, levels = c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"))
  jlims <- range(jsub$exprs)
  jmid <- min(jlims) + (max(jlims) - min(jlims)) / 2
  
  m.lines <- ggplot(trajs.sum, aes(x = pos, y = exprs)) + geom_line() + facet_grid(lambda.bin ~ mark) + 
    theme_bw() + 
    theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(paste("Trajectory:", jtraj, jstr))
  m1 <- ggplot(jsub, aes(x = pos / 1e6, y = reorder(lambda.bin, desc(lambda.bin)), fill = exprs)) + 
    geom_tile() + 
    theme_bw() + theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
    scale_fill_gradient(high = "darkblue", low = "gray85") + 
    facet_wrap(~mark, ncol = 1) + 
    ggtitle(paste(jtraj, jstr)) + 
    ylab("Trajectory") + 
    xlab("Position (MB)")
  print(m.lines)
  print(m1)

}
if (make.plots){
  dev.off()
}

if (make.plots){
  pdf(paste0("~/data/scchic/pdfs/variance_over_pseudotime_plots_primetime.Figure4.TrajsSD.", Sys.Date(), ".pdf"), useDingbats = FALSE)
}
  # Prepare dat -------------------------------------------------------------

  # get groupings by trajecotry
  
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
  
  # label trajectory and lambda
  
  # add lambda
  trajs.long <- lapply(jtrajs, function(jtraj){
    trajs.tmp <- lapply(trajs.mixed, function(x) x[[jtraj]]) %>%
      bind_rows() %>%
      rowwise() %>%
      mutate(mark = strsplit(cell, "_")[[1]][[2]]) %>%
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
                       axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom") +  
    scale_color_identity() + 
    xlab("Pseudotime") + ylab("Genome-wide SD") 
  
  m.nofacet <- ggplot(jsub, aes(x = lambda, y = cell.sd, color = jcol, group = traj)) + 
    facet_wrap(~mark, nrow = 1) + 
    geom_point(alpha = 0.3) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                       axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom") +  
    scale_color_identity() + 
    xlab("Pseudotime") + ylab("Genome-wide SD") 
  
  m.lst <- lapply(jtrajs, function(jtraj){
    m <- ggplot(jsub.trajfilt, aes(x = lambda, y = cell.sd, color = jcol, group = traj)) + 
      facet_wrap(~mark, nrow = 1) + 
      geom_point(alpha = 0.3) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                         axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom") + 
      ggtitle(jtraj) + 
      scale_color_identity() + 
      xlab("Pseudotime") + ylab("Genome-wide SD") 
    return(m)
  })
  
  print(m.facet)
  print(m.nofacet)
  print(m.lst)

  if (make.plots){
    dev.off()
  }


