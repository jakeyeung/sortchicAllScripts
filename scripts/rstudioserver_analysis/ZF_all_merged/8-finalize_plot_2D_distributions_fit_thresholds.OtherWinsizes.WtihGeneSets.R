# Jake Yeung
# Date of Creation: 2020-05-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/8-finalize_plot_2D_distributions_fit_thresholds.OtherWinsizes.WtihGeneSets.R
# description

rm(list=ls())

library(ggrastr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(preprocessCore)

library(mixtools)

library(JFuncs)


AssignQuadrant <- function(active.s2n.fc, repress.s2n.fc){
  if (any(is.na(c(active.s2n.fc, repress.s2n.fc)))){
    return(NA)
  }
  if (active.s2n.fc > 0 & repress.s2n.fc < 0){
    jquad <- "quad4_activeHigh_repressLow"
  } else if (active.s2n.fc <= 0 & repress.s2n.fc <= 0) {
    jquad <- "quad3_activeLow_repressLow"
  } else if (active.s2n.fc < 0 & repress.s2n.fc > 0){
    jquad <- "quad2_activeLow_repressHigh"
  } else if (active.s2n.fc >= 0 & repress.s2n.fc >= 0){
    jquad <- "quad1_activeHigh_repressHigh"
  } else {
    print("Cannot assign:")
    print(paste(active.s2n.fc, repress.s2n.fc))
  }
  return(jquad)
}




FitThresholdMM <- function(x, jlambda = c(0.5, 0.5), jmu = c(2, 8), jk = 2, jsigma = c(2, 2), jthres = 0.5, show.plot = TRUE, plot.title = ""){
  jout <- mixtools::normalmixEM(x, lambda = jlambda, mu = jmu, k = jk, sigma = jsigma)
  (xline <- min(jout$x[which(jout$posterior[, 1] < jthres)]))
  if (show.plot){
    plot.mixEM(jout, whichplots = 2, main2 = plot.title); abline(v = xline, lwd=3, lty=2, col = 'blue')
  }
  return(list(mm = jout, xline = xline))
}

jalpha <- 0.15
jdotsize <- 5

# try the hardest one first, K27me3 sqrt for Bcells?
jlambda <- c(0.5, 0.5)
jmu <- c(1, 3)
jk <- 2
jsigma <- c(1, 1)
jthres <- 0.5
mincounts <- 3

jwinsize <- 20000L

jpos <- "bottom"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/pseudobulk_analysis.integrated_analysis.TSS.plot_2D_distributions.OtherWinsizes"
dir.create(outdir)
outname <- paste0("WKM_2D_active_vs_repressed_fit_thresholds.", Sys.Date(), ".winsize_", jwinsize, ".WithGeneSets.pdf")

outf <- file.path(outdir, outname)

make.plots <- TRUE
if (make.plots){
  pdf(outf, useDingbats = FALSE)
}


# Init inputs -------------------------------------------------------------

downsample.reads <- TRUE
downsample.cells <- FALSE
# jdate <- "2020-05-04"  # downsample.cells FALSE did nothing
jdate <- "2020-05-09"  # downsample.cells FALSE now does not down sample 
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/from_rstudio/pseudobulk_analysis_all/pseudobulk_analysis.integrated_analysis.TSS.readsDS.likeBM.OtherWinsizes"

# inprefix <- "integrated_analysis.2020-05-09.UseTSSfromH3K4me3.DSreads_TRUE.DScells_FALSE.likeBM.winsize_10000.RData"
inprefix <- paste0("integrated_analysis.", jdate, ".UseTSSfromH3K4me3.DSreads_", downsample.reads, ".DScells_", downsample.cells, ".likeBM.winsize_", jwinsize)
inrdata <- paste0(inprefix, ".RData")
infrdata <- file.path(indir, inrdata)

assertthat::assert_that(file.exists(infrdata))

load(infrdata, v=T)


# Rename objects? ---------------------------------------------------------

jmerged <- jmerged.annot


# ggrastr::ggrastr::geom_point_rast_rast()

# Plot the active vs K27me3 across pbulks  --------------------------------

boxplot(jmerged %>% dplyr::select(c(H3K4me1, H3K4me3, H3K27me3)))

ggplot(jmerged, aes(x = sqrt(H3K4me3), y = sqrt(H3K27me3))) + 
  ggrastr::geom_point_rast(size = jdotsize, alpha = jalpha) + 
  facet_wrap(~cluster) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = jpos) + 
  geom_density2d()

# show the density plots find the background levels 

jpos <- "bottom"
# sqrt, log, or linear?
m.sqrt <- lapply(jmarks, function(jmark){
  print(jmark)
  m <- ggplot(jmerged, aes_string(x = paste0("sqrt(", jmark, ")"), fill = "cluster")) + geom_density(alpha=0.3) + 
    # facet_wrap(~cluster, ncol = 1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = jpos) + 
    ggtitle(jmark)
})

m.log <- lapply(jmarks, function(jmark){
  print(jmark)
  m <- ggplot(jmerged, aes_string(x = paste0("log2(", jmark, " + 1)"), fill = "cluster")) + geom_density(alpha=0.3) + 
    # facet_wrap(~cluster, ncol = 1) +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = jpos) + 
    ggtitle(jmark)
})


m.linear <- lapply(jmarks, function(jmark){
  print(jmark)
  m <- ggplot(jmerged, aes_string(x = jmark, fill = "cluster")) + geom_density(alpha=0.3) + 
    # facet_wrap(~cluster, ncol = 1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = jpos) + 
    ggtitle(jmark)
})

m.fracs <- lapply(jmarks, function(jmark){
  print(jmark)
  m <- ggplot(jmerged, aes_string(x = paste0(jmark, "/ sum(", jmark, ")"), fill = "cluster")) + geom_density(alpha=0.3) + 
    # facet_wrap(~cluster, ncol = 1) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = jpos) + 
    ggtitle(jmark)
})


for (jmark in jmarks){
  multiplot(m.linear[[jmark]], m.sqrt[[jmark]] + geom_vline(xintercept = 7), m.log[[jmark]], cols = 3)
}


# Find background level by fitting Gaussian mixture model  ----------------



# x <- sqrt(subset(jmerged, cluster == "Bcells")$H3K27me3)
# 
# mm.out <- FitThresholdMM(x, plot.title="Bcells,H3K27me3")
# 
jclsts <- unique(jmerged$cluster)
names(jclsts) <- jclsts

# plot all using sqrt

set.seed(0)

dat.mm <- lapply(jmarks, function(jmark){
  jdat.lst <- lapply(jclsts, function(jclst){
    print(paste(jmark, jclst))
    jsub <- subset(jlong.lst[[jmark]], cluster == jclst)
    xcounts <- jsub$counts
    xcounts <- xcounts[which(xcounts > mincounts)]
    jout <- FitThresholdMM(x = sqrt(xcounts), plot.title = paste(jmark, jclst))
    # create data frame?
    jdat <- data.frame(thres.sqrt = jout$xline, mark = jmark, cluster = jclst, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  mutate(thres = thres.sqrt ^ 2,
         thres.log2 = log2(thres + 1)) 

# plot distributions with threshlds

m.sqrt.thres <- lapply(jmarks, function(jmark){
  print(jmark)
  jsub.dat <- jlong.lst[[jmark]]
  jthres.dat <- subset(dat.mm, mark == jmark)
  m <- ggplot(jsub.dat, aes(x = sqrt(counts), fill = cluster)) + 
    geom_density(alpha=0.3) +  
    geom_vline(mapping = aes(xintercept = sqrt(thres), color = cluster), data = jthres.dat) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") +  
    facet_wrap(~cluster, ncol = 1) + 
    ggtitle(jmark)
})
print(m.sqrt.thres)

# calculate signal to noise ratio 
jlong.thres.lst <- lapply(jmarks, function(jmark){
  jsub <- jlong.lst[[jmark]]
  jsub$mark <- jmark
  # calculate signal to noise ratio
  jsub.merge <- left_join(jsub, subset(dat.mm, mark == jmark)) %>%
    rowwise() %>%
    mutate(s2n.fc = log2(counts + 1) - thres.log2,
           s2n.diff = counts - thres,
           s2n.diffsqrt = sqrt(counts) - thres.sqrt)
})

# plot s2n for one mark across celltypes
for (jmark in jmarks){
  print(jmark)
  m <- ggplot(jlong.thres.lst[[jmark]], aes(x = s2n.diffsqrt, fill = cluster)) + geom_density(alpha=0.3) +
    
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
    ggtitle(jmark) + geom_vline(xintercept = 0, linetype = 'dotted') + 
    facet_wrap(~cluster, ncol = 1)
  print(m)
}

# plot K4me3 vs K27me3 distribution?

jmerged.s2n.fc <- reshape2::dcast(as.data.frame(jlong.thres.lst %>% bind_rows()), formula = "bin + cluster ~ mark", value.var = "s2n.fc")
jmerged.s2n.sqrt <- reshape2::dcast(as.data.frame(jlong.thres.lst %>% bind_rows()), formula = "bin + cluster ~ mark", value.var = "s2n.diffsqrt")
jmerged.s2n.lin <- reshape2::dcast(as.data.frame(jlong.thres.lst %>% bind_rows()), formula = "bin + cluster ~ mark", value.var = "s2n.diff")

ggplot(jmerged.s2n.fc, aes(x = H3K4me1, y = H3K27me3, color = cluster)) + 
  ggrastr::geom_point_rast(size = jdotsize, alpha = jalpha) + 
  facet_wrap(~cluster, ncol = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  geom_density2d(color = "black") + 
  ggtitle("Signal2Noise log2FoldChange")

ggplot(jmerged.s2n.fc, aes(x = H3K4me3, y = H3K27me3, color = cluster)) + 
  ggrastr::geom_point_rast(size = jdotsize, alpha = jalpha) + 
  facet_wrap(~cluster, ncol = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  geom_density2d(color = "black") + 
  ggtitle("Signal2Noise log2FoldChange")


ggplot(jmerged.s2n.sqrt, aes(x = H3K4me1, y = H3K27me3, color = cluster)) + 
  ggrastr::geom_point_rast(size = jdotsize, alpha = jalpha) + 
  facet_wrap(~cluster, ncol = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  geom_density2d(color = "black") + 
  ggtitle("Signal2Noise SquareRoot")

ggplot(jmerged.s2n.sqrt, aes(x = H3K4me3, y = H3K27me3, color = cluster)) + 
  ggrastr::geom_point_rast(size = jdotsize, alpha = jalpha) + 
  facet_wrap(~cluster, ncol = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  geom_density2d(color = "black") + 
  ggtitle("Signal2Noise SquareRoot Diff")

ggplot(jmerged.s2n.sqrt, aes(x = H3K4me3, y = H3K27me3, color = cluster)) + 
  ggrastr::geom_point_rast(size = jdotsize, alpha = jalpha) + 
  facet_wrap(~cluster, ncol = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  geom_density2d(color = "black") + 
  ggtitle("Signal2Noise SquareRoot Diff")

ggplot(jmerged.s2n.lin, aes(x = H3K4me1, y = H3K27me3, color = cluster)) + 
  ggrastr::geom_point_rast(size = jdotsize, alpha = jalpha) + 
  facet_wrap(~cluster, ncol = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  geom_density2d(color = "black") + 
  ggtitle("Signal2Noise Linear Diff")

ggplot(jmerged.s2n.lin, aes(x = H3K4me3, y = H3K27me3, color = cluster)) + 
  ggrastr::geom_point_rast(size = jdotsize, alpha = jalpha) + 
  facet_wrap(~cluster, ncol = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
  geom_density2d(color = "black") + 
  ggtitle("Signal2Noise Linear Diff")

# Plot as boxplots maybe cleanest  ----------------------------------------

ggplot(jlong.thres.lst %>% bind_rows(), aes(x = cluster, y = s2n.diffsqrt, fill = mark)) +
  geom_violin() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_hline(yintercept = 0, linetype = "dotted")


# Count fractions in the four quadrants?  ---------------------------------


jcountgenes.k4me3 <- jmerged.s2n.fc %>% 
  bind_rows() %>%
  rowwise() %>%
  mutate(quadrant = AssignQuadrant(H3K4me3, H3K27me3)) %>%
  filter(!is.na(quadrant)) %>%
  group_by(cluster, quadrant) %>%
  summarise(ngenes = length(bin)) %>%
  group_by(cluster) %>%
  mutate(fracgenes = ngenes / sum(ngenes))
ggplot(jcountgenes.k4me3, aes(x = cluster, y = fracgenes, group = quadrant, fill = quadrant)) + 
  geom_col(position = "dodge") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle("HSPC enriched in activeHigh (K4me3) and repressHigh (K27me3) relative to other celltypes")

jcountgenes.k4me1 <- jmerged.s2n.fc %>% 
  bind_rows() %>%
  rowwise() %>%
  mutate(quadrant = AssignQuadrant(H3K4me1, H3K27me3)) %>%
  filter(!is.na(quadrant)) %>%
  group_by(cluster, quadrant) %>%
  summarise(ngenes = length(bin)) %>%
  group_by(cluster) %>%
  mutate(fracgenes = ngenes / sum(ngenes))
ggplot(jcountgenes.k4me1, aes(x = cluster, y = fracgenes, group = quadrant, fill = quadrant)) + 
  geom_col(position = "dodge") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ggtitle("HSPC enriched in activeHigh (K4me1) and repressHigh (K27me3) relative to other celltypes")

# K4me1 is highest for HSPC high
jcounts.check1 <- jmerged.s2n.fc %>%
  bind_rows() %>%
  group_by(cluster) %>%
  summarise(ngenes = length(which(H3K4me1 > 0)))
print(jcounts.check1)

jcounts.check2 <- jmerged.s2n.fc %>%
  bind_rows() %>%
  group_by(cluster) %>%
  summarise(ngenes = length(which(H3K27me3 > 0)))
print(jcounts.check2)


# Save outputs ------------------------------------------------------------

if (file.exists(outf.rdata)){
  print(paste("Skipping saving rdata, exists:", outf.rdata))
} else {
  save(jlong.thres.lst, mm.out, file = outf.rdata)
}



if (make.plots){
  dev.off()
}


