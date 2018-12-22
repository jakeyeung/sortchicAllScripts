# Jake Yeung
# plot_peak_counts_for_thresholds.R
# Make exploratory figures for figuring out proper thresholds for excluding cells and peaks
# 2018-12-22

library(dplyr)
library(ggplot2)

args <- commandArgs(trailingOnly=TRUE)

source("scripts/Rfunctions/GetMetaData.R")

inf <- args[[1]]  # count matrix object from Rsubreads
outdir <- args[[2]]

load(inf, v=T)
count.mat <- count.dat$counts

outpath <- file.path(outdir, "count_matrix_plots.pdf")

# Plot output -------------------------------------------------------------

dat.meanvar <- data.frame(Sum = Matrix::rowSums(count.mat), 
                          Mean = Matrix::rowMeans(count.mat),
                          Var = apply(count.mat, 1, var),
                          peak = rownames(count.mat), 
                          stringsAsFactors = FALSE)
dat.meanvar <- dat.meanvar %>%
  rowwise() %>%
  mutate(CV = sqrt(Var) / Mean,
         peaksize = GetPeakSize(peak))

pdf(outpath, useDingbats=FALSE)

  # maybe a few suspicious peaks?
  ggplot(dat.meanvar, aes(x = peaksize, y = Sum)) + geom_point(alpha = 0.1) +
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_y_log10()
  
  ggplot(dat.meanvar, aes(x = peaksize)) + geom_density() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  ggplot(dat.meanvar, aes(x = log2(Mean), y = log2(CV), size = peaksize)) + geom_point(alpha = 0.25) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_abline(slope = -0.5) + ggtitle("Includes some suspicious peaks")
  
  # plot cell sizes 
  plot(density(colSums(count.mat)))
  plot(hist(colSums(count.mat), breaks = 75))

dev.off()

