# Jake Yeung
# Date of Creation: 2019-05-03
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/variance_over_pseudotime_heatmaps_zoomed_in.R
# Zoom in on genomic regions 


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

library(GGally)

source("scripts/Rfunctions/VariabilityFunctions.R")
source("scripts/Rfunctions/Aux.R")
source("scripts/Rfunctions/TrajFunctions.R")





# Make as function --------------------------------------------------------


jtraj <- "granu"


# Load data  --------------------------------------------------------------

# inf.trajs <- "/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-09.RData"
inf.trajs <- "/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-11.RData"
inf.trajs <- "/Users/yeung/data/scchic/robjs/trajectory_from_spring_2019-04-15.RData"
assertthat::assert_that(file.exists(inf.trajs))
load(inf.trajs, v=T)

head(trajs.spring[[4]][[1]])

inf.dat <- "/Users/yeung/data/scchic/robjs/TFactivity_genelevels_objects_build95.allmarks_reorient_WithTrajs.WithColnamesLst.2019-04-04.RData"
assertthat::assert_that(file.exists(inf.dat))
load(inf.dat, v=T)


trajs.mixed.lst <- GetTrajMixed()
trajs.spring <- trajs.mixed.lst$trajs.mixed
trajs.mixed <- trajs.mixed.lst$trajs.mixed

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

trajs.sum <- lapply(jstrs, function(jstr) GetTrajSum(tm.result.lst = tm.result.lst, trajs.mixed = trajs.mixed, jmarks, jstr = jstr, jfac = jfac, jpseudo = jpseudo, jtraj = jtraj)) %>%
  bind_rows()

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

outf <- "~/data/scchic/robjs/trajs_sum_variance.rds"
if (!file.exists(outf)){
  saveRDS(trajs.sum, file = outf)
}

pdf(file = paste0("~/data/scchic/pdfs/variance_over_pseudotime_slopes_and_correlations_genomewide.", Sys.Date(), ".pdf"), useDingbats = FALSE)
  # show mutual exclusiivity between the two repressive marks 
  for (jchromo in jstrs){
    m.slopes <- ggplot(jsub.fits %>% filter(chromo == jchromo), aes(x = pos / 10^6, y = reorder(mark, desc(mark)), fill = slope)) + geom_tile() + 
      theme_bw() + theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
      scale_fill_gradient2(low = "darkred", mid = "gray85", high = "darkblue", midpoint = 0, name = "Slope [A.U.]") + 
      ylab("") + 
      xlab("Position (MB)") + 
      facet_wrap(~chromo)
    print(m.slopes)
  }
  # plot mutual exclusivity
  m.cor <- ggpairs(jsub.wide %>% select(-pos, -chromo),
                   lower = list(continuous = HexPlot)) + theme_classic() + ggtitle(jstr)
  print(m.cor)
dev.off()


