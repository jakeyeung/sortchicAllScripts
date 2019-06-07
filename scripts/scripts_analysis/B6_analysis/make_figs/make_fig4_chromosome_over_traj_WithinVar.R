# Jake Yeung
# Date of Creation: 2019-06-03
# File: ~/projects/scchic/scripts/scripts_analysis/B6_analysis/make_figs/make_fig4_chromosome_over_traj_WithinVar.R
# Use WithinVar as signal
# find colors https://menugget.blogspot.com/2011/11/define-color-steps-for-colorramppalette.html#more


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


inf.var <- paste0("/Users/yeung/data/scchic/robjs/B6_objs/variance_decomposed_within_across.RData")
assertthat::assert_that(file.exists(inf.var))
inf.traj <- paste0("/Users/yeung/data/scchic/robjs/B6_objs/traj_objs_all_marks.Rdata")
assertthat::assert_that(file.exists(inf.traj))
indir <- "/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6"
assertthat::assert_that(dir.exists(indir))

pdfout <- paste0("/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory/variance_over_pseudotime_fewer_trajs.withinvar.", Sys.Date(), ".pdf")

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


# Load variance -----------------------------------------------------------

load(inf.var, v=T)
chromo.constant.sum <- chromo.constant %>%
  summarise(nbins = sum(nbins))

# Load data  --------------------------------------------------------------


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
# jmark <- jmarks[["H3K4me1"]]
infs <- lapply(jmarks, function(jmark) list.files(indir, pattern = paste0(jmark, ".RData"), full.names = TRUE))
tm.result.lst <- lapply(infs, LoadGetTmResult)


# Load trajectories -------------------------------------------------------


load(inf.traj, v=T)



# Plot variance over pseudotime  ------------------------------------------



pdf(pdfout, useDingbats = FALSE)

# pdf("/Users/yeung/data/scchic/pdfs/B6_figures/variance_over_trajectory/variance_over_pseudotime_fewer_trajs.withinvar.2019-06-03.debug.pdf", useDingbats=FALSE)

# variance over pseudotime 
colhash <- GetTrajColors(as.hash = TRUE, add.mega = FALSE)
jsub <- cells.var.merged2 %>% filter(!is.na(traj))
jsub$jcol <- sapply(jsub$traj, function(x) colhash[[x]])
m.var1 <- ggplot(jsub %>% filter(varname == "cell.var.within.sum"), aes(x = lambda, y = varval / chromo.constant.sum$nbins, color = jcol)) + geom_point(alpha = 0.2) + facet_wrap(~mark, nrow = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_continuous(breaks = c(0, 1)) + 
  xlab("Pseudotime") + 
  ylab(expression('Variance ' * group("[", log[2] * ("signal")^{2}, "]"))) + 
  scale_color_identity()
print(m.var1)
# dev.off()

m.var.smooth <- ggplot(jsub %>% filter(varname == "cell.var.within.sum"), aes(x = lambda, y = varval / chromo.constant.sum$nbins, color = jcol)) + geom_point(alpha = 0.2) + facet_wrap(~mark, nrow = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_smooth(method = "lm", se = FALSE, size = 2) + 
  scale_x_continuous(breaks = c(0, 1)) + 
  xlab("Pseudotime") +
  ylab(expression('Variance ' * group("[", log[2] * ("signal")^{2}, "]"))) + 
  scale_color_identity()
print(m.var.smooth)

# umap
dat.merged <- left_join(bind_rows(dat.umap.long.trajs), cells.var.merged2 %>% filter(varname == "cell.var.within.sum") %>% dplyr::select(cell, varval)) %>%
  mutate(cell.within.var = varval / chromo.constant.sum$nbins)

# library("colorspace")
# pal <- colorspace::choose_palette()

jtitle <- ""
jcname <- "cell.within.var"
xvar <- "umap1"
yvar <- "umap2"
# steps <- c("blue4", "cyan", "white", "yellow", "red4")
steps <- c("gray85", "gray50", "gray40",  "cyan")
pal <- color.palette(steps, c(6, 2, 40), space="rgb")
jsize <- 4

for (jmark in jmarks){
  # m1 <- PlotXYWithColor(dat.merged %>% filter(mark == jmark), xvar = "umap1", yvar = "umap2", cname = "cell.within.var", jcol.low = "gray90", jcol.mid = "gray50", jcol = "cyan", jsize = 6)
  # print(m1)
  dat.merged$neg.cell.within.var <- -1 * dat.merged$cell.within.var
  dat.merged <- RankOrder(dat.merged, cname = "neg.cell.within.var", out.cname = "orderrank")
  m1 <- ggplot(dat.merged %>% filter(mark == jmark), aes_string(x = xvar, y = yvar, col = jcname, order = "orderrank")) +
    ggrastr::geom_point_rast(size = jsize) +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom",
                       axis.ticks=element_blank(),
                       axis.text.x=element_blank(),
                       axis.text.y=element_blank(),
                       panel.border=element_blank())  +
    xlab("") + ylab("") + ggtitle(jmark) +
    # scale_color_gradientn(colours = RColorBrewer::brewer.pal(20, "YlGnBu"))
    scale_color_gradientn(colours = pal(50))
    # scale_color_gradient2(low = "gray85", high = "cyan", mid = "gray50", breaks = c(0.5, 2, 4), midpoint = 2)
  print(m1)
}

for (jtrajname in jtrajs){
  m.var <- ggplot(cells.var.merged2 %>% filter(traj == jtrajname), aes(x = lambda, y = varval / chromo.constant.sum$nbins, color = varname, group = varname)) + geom_point(alpha = 0.2) + facet_wrap(~mark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_smooth(method = "lm", se = FALSE, size = 2) + 
    xlab("Pseudotime") + ylab("Variance (log2 signal ^ 2)") + 
    scale_color_discrete(name = "", labels = c("VarGlobal", "VarAcross", "VarWithin")) + 
    ggtitle(jtrajname)
  print(m.var)
  
  ssdiff <- sweep(x = cells.chromo.mean.wide.mat, MARGIN = 2, STATS = dat.mat.global.mean, FUN = "-") ^ 2
  ssdiff <- sweep(x = ssdiff, MARGIN = 1, STATS = chromo.constant$nbins, FUN = "*")
  chromo.mean.dat <- data.frame(chromo.var = unlist(data.frame(ssdiff)), 
                                chromo = rep(rownames(ssdiff), ncol(ssdiff)), 
                                cell = rep(colnames(ssdiff), each = nrow(ssdiff)), stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(mark = strsplit(cell, "-")[[1]][[4]])
  m.withinbychromo <- ggplot(cells.var.chromo.within, aes(x = log10(cell.var.within), fill = mark)) + geom_density(alpha = 0.25) + facet_wrap(~label)
  m.acrossvar <- ggplot(chromo.mean.dat, aes(x = log10(chromo.var), fill = mark)) + geom_density(alpha = 0.25) + facet_wrap(~chromo) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m.withinbychromo)
  print(m.acrossvar)
}
  dev.off()


# Plot umap with variance  ------------------------------------------------


