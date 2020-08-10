# Jake Yeung
# Date of Creation: 2020-06-21
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/downstream_poisson_model_analysis_show_2D_cloud.with_CI.OtherNorms.R
# description


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(ggrepel)
library(scchicFuncs)
library(JFuncs)

make.plots <- TRUE

# jnorm <- "ncuts.inbins"
jnorm <- "ncuts.alltss"
bsize <- 10000

jdate <- "2020-06-06"
jdate2 <- "2020-06-08"

fewer.k27me3 <- TRUE
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.setup_for_poisson_regression"
jprefix <- file.path(indir, paste0("integrated_analysis_3_marks.stringentDE.WithHighLowExprs.singlecells.fewerk27me3_", fewer.k27me3, ".forPoissonRegression.CountR1only.", jdate))

# outfits <- file.path(indir, paste0("fit_poisson_model_on_TSS.", jdate, ".RData"))
fitsdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits.OtherNorm.DiffDists"
# jprefix <- paste0("fit_poisson_model_on_TSS.MouseBM.NormMeth_", jnorm, ".bsize_", bsize)
jprefix <- paste0("fit_poisson_model_on_TSS.MouseBM.NormMeth_", jnorm, ".bsize_", bsize, ".CleanUpErythss")
fitprefix <- file.path(fitsdir, jprefix)
infit <- paste0(fitprefix, ".RData")
infit.wrangled <- paste0(fitprefix, ".DownstreamWrangled.RData")
# infit.ci <- file.path(fitsdir, paste0("fit_poisson_model_on_TSS.MouseBM.NormMeth_", jnorm, ".CI.bsize_", bsize, ".DownstreamWrangled.RData"))
infit.ci <- file.path(fitsdir, paste0("fit_poisson_model_on_TSS.MouseBM.NormMeth_", jnorm, ".CI.bsize_", bsize, ".CleanUpErythss.DownstreamWrangled.RData"))

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/integrated_analysis_poisson_and_2D_clouds"
assertthat::assert_that(dir.exists(outdir))

# pdfout <- file.path(outdir, paste0(jprefix, ".Downstream2DClouds.", Sys.Date(), ".WithCI.pdf"))
# mixedbinsout <- file.path(outdir, paste0(jprefix, ".Downstream2DClouds.", Sys.Date(), ".WithCI.MixedBins.txt"))

pdfout <- file.path(outdir, paste0(jprefix, ".Downstream2DClouds.", Sys.Date(), ".WithCI.CleanUpErythss.pdf"))
mixedbinsout <- file.path(outdir, paste0(jprefix, ".Downstream2DClouds.", Sys.Date(), ".WithCI.MixedBins.CleanUpErythss.txt"))

# assertthat::assert_that(!file.exists(pdfout))
# assertthat::assert_that(!file.exists(mixedbinsout))

load(infit.wrangled, v=T)
load(infit.ci, v=T)

jfits.long <- jfits.long %>%
  group_by(bin) %>%
  filter(abs(logLambda) < 10) %>%
  filter(abs(logintercept) < 20)

# add arrows? check for Bcell-specific genes vs Eryth+Granu specific genes
print(unique(fits.bygenesets.long$geneset))
print(unique(fits.bygenesets.long$cluster))

ggplot(jfits.long, aes(x = logLambda, fill = cluster)) + geom_density(alpha = 0.25) + facet_wrap(~mark) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = 0, linetype = "dotted")


# Make into matrix and plot 2D matrix  ------------------------------------


# Plot the 2D cloud -------------------------------------------------------


ctypes.end <- list("ClusterGranulocytes", "ClusterBcells", "ClusterErythroblasts"); names(ctypes.end) <- ctypes.end
gset.specs <- list("Neutrophil", "Bcell", "Erythroblast", "HighExprs"); names(gset.specs) <- ctypes.end
gset.others <- list(c("Bcell", "Erythroblast"), c("Erythroblast", "Neutrophil"), c("Bcell", "Neutrophil"), c("LowExprs")); names(gset.others) <- ctypes.end
gsets.differentiated <- c("Neutrophil", "Bcell", "Erythroblast")
gset.hsc <- "HSCs"

ctype.end <- ctypes.end[[1]] 

jscale2 <- 0.2
low.cutoff <- -11
# take top 5% 
jprob <- 0.8

# fc.max <- 20
if (make.plots){
  # pdf(pdfout, width = 1020/72, height = 815/72, useDingbats = FALSE)
  pdf(pdfout, useDingbats = FALSE)
}

bins.in.gset <- unique(subset(fits.bygenesets.long, geneset %in% gsets.differentiated)$bin)

# plot genomewide 

jfits.long.tmp <- jfits.long %>%
  ungroup() %>%
  mutate(mark = factor(mark, levels = c("H3K27me3", "H3K4me1", "H3K4me3")),
         cluster = factor(as.character(cluster), levels = c("ClusterErythroblasts", "ClusterGranulocytes", "ClusterBcells")))

ggplot(jfits.long.tmp, aes(x = logLambda, fill = mark)) + geom_density() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark) + 
  geom_vline(xintercept = 0, linetype = "dotted")

ggplot(jfits.long.tmp, aes(x = logLambda, fill = cluster)) + geom_density(alpha = 0.33) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark) + 
  geom_vline(xintercept = 0, linetype = "dotted")

ggplot(jfits.long.tmp, aes(x = logLambda, fill = cluster)) + geom_density(alpha = 0.33) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_grid(cluster~mark) + 
  geom_vline(xintercept = 0, linetype = "dotted")



# plot against all genesets 
cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

fits.bygenesets.long.tmp <- fits.bygenesets.long %>%
  ungroup() %>%
  mutate(cluster = factor(cluster, levels = c("ClusterErythroblasts", "ClusterGranulocytes", "ClusterBcells")))
for (jmark in jmarks){
  m.gset <- ggplot(fits.bygenesets.long.tmp %>% filter(abs(logLambda) < 5 & mark == jmark), aes(x = logLambda, fill = geneset)) + geom_density(alpha = 0.5) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_fill_manual(values = cbPalette) + 
    facet_grid(geneset ~ cluster) + geom_vline(xintercept = 0, linetype = "dotted") + 
    ggtitle(jmark)
  print(m.gset)
}
  

lapply(ctypes.end, function(ctype.end){
 
   
  print(ctype.end)
  
  gset.spec <- gset.specs[[ctype.end]]
  gset.other <- gset.others[[ctype.end]]
  
  
  jtitle <- paste0("HSCs -> ", ctype.end)

  jfits.mat.ints <- reshape2::dcast(data = jfits.long %>% filter(cluster == ctype.end), 
                               formula = bin + gene + ens ~ mark, value.var = "logintercept")
  jfits.mat.logfcs <- reshape2::dcast(data = jfits.long %>% filter(cluster == ctype.end), 
                               formula = bin + gene + ens ~ mark, value.var = "logLambda") %>%
    dplyr::rename(H3K4me1.fc = H3K4me1,
                  H3K4me3.fc = H3K4me3,
                  H3K27me3.fc = H3K27me3)
  
  
  jfcs.all <- jfits.mat.logfcs %>%
    rowwise() %>%
    left_join(., subset(jfits.mat.ints, select = c(bin, H3K4me1, H3K4me3, H3K27me3)), by = "bin")
  jfcs.mixed <- subset(jfcs.all, H3K4me3 >= quantile(H3K4me3, probs = jprob, na.rm = TRUE) & H3K27me3 >= quantile(H3K27me3, prob = jprob, na.rm = TRUE))
  nbins.mixed <- nrow(jfcs.mixed)
  
  
  # plot mixed 
  jfcs.all.mixed_vs_all <- jfcs.all %>% filter(H3K4me3 >= low.cutoff & H3K27me3 >= low.cutoff) %>%
    rowwise() %>%
    mutate(gset = ifelse(bin %in% jfcs.mixed$bin, "Mixed", "zNotMixed"))
  
  m.cloud.arrows <- ggplot(jfcs.all.mixed_vs_all, aes(x = H3K4me3, y = H3K27me3, color = gset)) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_segment(mapping = aes(xend = H3K4me3 + jscale2 * H3K4me3.fc, yend = H3K27me3 + jscale2 * H3K27me3.fc),
                 arrow = arrow(length=unit(0.25,"cm"), ends = "last"), alpha = 0.8, size = 0.1) + 
    ggtitle(jtitle, paste("genes levels greater than", low.cutoff))
  print(m.cloud.arrows)
  
  xrange.cloud.log <- ggplot_build(m.cloud.arrows)$layout$panel_scales_x[[1]]$range$range
  yrange.cloud.log <- ggplot_build(m.cloud.arrows)$layout$panel_scales_y[[1]]$range$range
  
  m.cloud.arrows.mixed <- ggplot(jfcs.all.mixed_vs_all %>% filter(gset == "Mixed"), aes(x = H3K4me3, y = H3K27me3, color = gset)) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_segment(mapping = aes(xend = H3K4me3 + jscale2 * H3K4me3.fc, yend = H3K27me3 + jscale2 * H3K27me3.fc),
                 arrow = arrow(length=unit(0.25,"cm"), ends = "last"), alpha = 0.8, size = 0.1) + 
    ggtitle(jtitle, paste("mixed states log, N:", nbins.mixed)) + 
    coord_cartesian(xlim = xrange.cloud.log, ylim = yrange.cloud.log)
  print(m.cloud.arrows.mixed)
  
  m.cloud.arrows.linear <- ggplot(jfcs.all.mixed_vs_all, aes(x = exp(H3K4me3), y = exp(H3K27me3), color = gset)) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_segment(mapping = aes(xend = exp(H3K4me3 + H3K4me3.fc), yend = exp(H3K27me3 + H3K27me3.fc)),
                 arrow = arrow(length=unit(0.25,"cm"), ends = "last"), alpha = 0.4, size = 0.1) + 
    ggtitle(jtitle, paste("mixed states linear, N:", nbins.mixed))
  print(m.cloud.arrows.linear)
  
  xrange.cloud.linear <- ggplot_build(m.cloud.arrows.linear)$layout$panel_scales_x[[1]]$range$range
  yrange.cloud.linear <- ggplot_build(m.cloud.arrows.linear)$layout$panel_scales_y[[1]]$range$range
  
  m.cloud.arrows.mixed.linear <- ggplot(jfcs.all.mixed_vs_all %>% filter(gset == "Mixed"), aes(x = exp(H3K4me3), y = exp(H3K27me3), color = gset)) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_segment(mapping = aes(xend = exp(H3K4me3 + H3K4me3.fc), yend = exp(H3K27me3 + H3K27me3.fc)),
                 arrow = arrow(length=unit(0.25,"cm"), ends = "last"), alpha = 0.4, size = 0.1) + 
    ggtitle(jtitle, paste("mixed states linear, N:", nbins.mixed)) + 
    coord_cartesian(xlim = xrange.cloud.linear, ylim = yrange.cloud.linear)
  print(m.cloud.arrows.mixed.linear)
  
  # show h3k4me3 vs h3k27me3
  m.fc.k4me3_vs_k27me3.mixed <- ggplot(jfcs.mixed, aes(x = H3K4me3.fc, y = H3K27me3.fc)) + 
    geom_point(alpha = 0.4, color = 'red') + 
    geom_density_2d(alpha = 0.8, color = "black") + 
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle, paste("mixed states, N:", nbins.mixed))
  print(m.fc.k4me3_vs_k27me3.mixed)
  
  xrange.fc <- ggplot_build(m.fc.k4me3_vs_k27me3.mixed)$layout$panel_scales_x[[1]]$range$range
  yrange.fc <- ggplot_build(m.fc.k4me3_vs_k27me3.mixed)$layout$panel_scales_y[[1]]$range$range
  
  # add CI
  jfcs.mixed.merge <- left_join(jfcs.mixed, jmat.fc.ci.lst[[ctype.end]], by = "bin")
  nbins2.mixed <- nrow(jfcs.mixed.merge)
  
  m.fc.k4me3_vs_k27me3.ci.mixed <- ggplot(jfcs.mixed.merge %>% filter(!is.na(H3K4me3.fc) & !is.na(H3K27me3.fc)), aes(x = H3K4me3.fc, y = H3K27me3.fc)) + 
    geom_errorbar(mapping = aes(ymin = H3K27me3.fc.lower, ymax = H3K27me3.fc.upper), alpha = 0.1, width = 0) + 
    geom_errorbarh(mapping = aes(xmin = H3K4me3.fc.lower, xmax = H3K4me3.fc.upper), alpha = 0.1, height = 0) + 
    geom_point(alpha = 0.4, color = 'red') + 
    geom_density_2d(alpha = 0.5, color = "black") + 
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle, paste("Nbins:", nbins2.mixed)) + 
    coord_cartesian(xlim = xrange.fc, ylim = yrange.fc)
  print(m.fc.k4me3_vs_k27me3.ci.mixed)
  
  
  if (ctype.end == ctypes.end[[1]]){
    
    hsc.bins <- unique(subset(fits.bygenesets.long, geneset == gset.hsc)$bin)
    assertthat::assert_that(length(hsc.bins) > 0)
    # show expression of HSPC-specific genes
    jfits.mat.ints.hsc <- jfits.mat.ints %>%
      mutate(hspc.spec.gene = bin %in% hsc.bins) %>%
      ungroup() %>%
      # arrange(desc(hspc.spec.gene)) %>%
      arrange(hspc.spec.gene) %>%
      filter(H3K4me3 >= low.cutoff & H3K27me3 >= low.cutoff)
    
    m.hsc <- ggplot(jfits.mat.ints.hsc, aes(x = H3K4me3, y = H3K27me3, color = hspc.spec.gene)) + 
      geom_point(alpha = 0.4) + 
      # facet_wrap(~hspc.spec.gene) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle("HSCs lambda, label by HSPC-specific or not")
    print(m.hsc)
    
    m.hsc.lin <- ggplot(jfits.mat.ints.hsc, aes(x = exp(H3K4me3), y = exp(H3K27me3), color = hspc.spec.gene)) + 
      geom_point(alpha = 0.4) + 
      # facet_wrap(~hspc.spec.gene) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle("HSCs lambda, label by HSPC-specific or not")
    print(m.hsc.lin)
    
    
    
    m <- ggplot(jfits.mat.ints, aes(x = H3K4me3, y = H3K27me3)) + 
      geom_point(alpha = 0.1) + 
      geom_density_2d() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle("HSCs lambda")
    print(m)
    
    m <- ggplot(jfits.mat.ints, aes(x = exp(H3K4me3), y = exp(H3K27me3))) + 
      geom_point(alpha = 0.1) + 
      geom_density_2d() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle("HSCs lambda")
    print(m)
    
    m <- ggplot(jfits.mat.ints, aes(x = exp(H3K4me1))) + 
      geom_density() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle("HSC lambda")
    print(m)
    
    m <- ggplot(jfits.mat.ints, aes(x = exp(H3K4me3))) + 
      geom_density() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle("HSC lambda")
    print(m)
    
    m <- ggplot(jfits.mat.ints, aes(x = exp(H3K27me3))) + 
      geom_density() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle("HSC lambda")
    print(m)
    
    m <- ggplot(jfits.mat.logfcs, aes(x = H3K4me3.fc, y = H3K27me3.fc)) + 
      geom_point(alpha = 0.1) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      geom_hline(yintercept = 0, color = 'blue') + geom_vline(xintercept = 0, color = 'blue') + 
      ggtitle(jtitle, "log Fold changes")
    print(m)
    
    m <- ggplot(jfits.mat.logfcs, aes(x = H3K4me1.fc, y = H3K27me3.fc)) + 
      geom_point(alpha = 0.1) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      geom_hline(yintercept = 0, color = 'blue') + geom_vline(xintercept = 0, color = 'blue') + 
      ggtitle(jtitle, "log Fold changes")
    print(m)
    
    m <- ggplot(jfits.mat.logfcs, aes(x = H3K4me1.fc)) + 
      geom_density() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      geom_vline(xintercept = 0, color = 'blue') + 
      ggtitle(jtitle, "log Fold changes")
    print(m)
    
    m <- ggplot(jfits.mat.logfcs, aes(x = H3K4me3.fc)) + 
      geom_density() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      geom_vline(xintercept = 0, color = 'blue') + 
      ggtitle(jtitle, "log Fold changes")
    print(m)
    
    m <- ggplot(jfits.mat.logfcs, aes(x = H3K27me3.fc)) + 
      geom_density() + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      geom_vline(xintercept = 0, color = 'blue') + 
      ggtitle(jtitle, "log Fold changes")
    print(m)
    
    
    
    # show arrows genomewide
    
    
    
    # show HighExprs and LowExprs genes
    jbins.high <- subset(fits.bygenesets.long, geneset %in% "HighExprs")$bin
    jbins.low <- subset(fits.bygenesets.long, geneset %in% "LowExprs")$bin
    
    # remove overlaps
    jbins.hl <- intersect(jbins.high, jbins.low)
    nbins.hl <- length(jbins.hl)
    assertthat::assert_that(length(jbins.hl) == 0)
    
    jfcs.hl <- subset(jfcs.all, bin %in% c(jbins.high, jbins.low)) %>%
      rowwise() %>%
      mutate(gset = ifelse(bin %in% jbins.high, "HighExprs", "LowExprs"))
    
    m.hl.linear <- ggplot(jfcs.hl, aes(x = exp(H3K4me3), y = exp(H3K27me3), color = gset)) + 
      geom_point(alpha = 0.1) +
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      facet_wrap(~gset) + ggtitle("Low and High Exprs Genes in Linear Scale")
    print(m.hl.linear)
    
    m.hl.log <- ggplot(jfcs.hl, aes(x = H3K4me3, y = H3K27me3, color = gset)) + 
      geom_point(alpha = 0.1) +
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      facet_wrap(~gset) + ggtitle("Low and High Exprs Genes in Log Scale")
    print(m.hl.log)
    
    m.cloud.arrows.hl.log <- ggplot(jfcs.hl, aes(x = H3K4me3, y = H3K27me3, color = gset)) + 
      facet_wrap(~gset) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      geom_segment(mapping = aes(xend = H3K4me3 + jscale2 * H3K4me3.fc, yend = H3K27me3 + jscale2 * H3K27me3.fc),
                   arrow = arrow(length=unit(0.25,"cm"), ends = "last"), alpha = 0.4, size = 0.1) + 
      ggtitle(jtitle, paste("High and Low Exprs genes log, N:", nbins.hl))
    print(m.cloud.arrows.hl.log)
    
    m.cloud.arrows.hl.linear <- ggplot(jfcs.hl, aes(x = exp(H3K4me3), y = exp(H3K27me3), color = gset)) + 
      facet_wrap(~gset) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      geom_segment(mapping = aes(xend = exp(H3K4me3 + H3K4me3.fc), yend = exp(H3K27me3 + H3K27me3.fc)),
                   arrow = arrow(length=unit(0.25,"cm"), ends = "last"), alpha = 0.4, size = 0.1) + 
      ggtitle(jtitle, paste("High and Low Exprs Genes linear, N:", nbins.hl))
    print(m.cloud.arrows.hl.linear)
    
    # summarize fold changes with arrows
    jfcs.hl.merge <- left_join(jfcs.hl, jmat.fc.ci.lst[[ctype.end]], by = "bin")
    nbins2.hl <- nrow(jfcs.hl.merge)
    
    # show h3k4me3 vs h3k27me3
    m.fc.k4me3_vs_k27me3.hl <- ggplot(jfcs.hl, aes(x = H3K4me3.fc, y = H3K27me3.fc, color = gset)) +
      facet_wrap(~gset) + 
      geom_point(alpha = 0.4) + 
      geom_density_2d(alpha = 0.8, color = "black") + 
      geom_vline(xintercept = 0) + 
      geom_hline(yintercept = 0) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle(jtitle, paste("mixed states, N:", nbins.hl))
    print(m.fc.k4me3_vs_k27me3.hl)
    
    xrange.fc <- ggplot_build(m.fc.k4me3_vs_k27me3.hl)$layout$panel_scales_x[[1]]$range$range
    yrange.fc <- ggplot_build(m.fc.k4me3_vs_k27me3.hl)$layout$panel_scales_y[[1]]$range$range
    
    m.fc.k4me3_vs_k27me3.ci.hl <- ggplot(jfcs.hl.merge %>% filter(!is.na(H3K4me3.fc) & !is.na(H3K27me3.fc)), 
                                         aes(x = H3K4me3.fc, y = H3K27me3.fc, color = gset)) + 
      facet_wrap(~gset) + 
      geom_errorbar(mapping = aes(ymin = H3K27me3.fc.lower, ymax = H3K27me3.fc.upper), alpha = 0.1, width = 0) + 
      geom_errorbarh(mapping = aes(xmin = H3K4me3.fc.lower, xmax = H3K4me3.fc.upper), alpha = 0.1, height = 0) + 
      geom_point(alpha = 0.4) + 
      geom_density_2d(alpha = 0.5, color = "black") + 
      geom_vline(xintercept = 0) + 
      geom_hline(yintercept = 0) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle(jtitle, paste("Nbins:", nbins2.hl)) + 
      coord_cartesian(xlim = xrange.fc, ylim = yrange.fc)
    print(m.fc.k4me3_vs_k27me3.ci.hl)
    
    # write mixed state bins to output
    fwrite(x = subset(jfcs.mixed.merge, select = c(bin, gene, ens, H3K4me1, H3K4me3, H3K27me3)), file = mixedbinsout, sep = "\t")
    
    # show levels of differentiated gene sets compared to overall 
    
    
    jfcs.all.annot <- jfcs.all %>%
      rowwise() %>%
      mutate(gset = ifelse(bin %in% bins.in.gset, "WillBeActivated", "zOther")) %>%
      ungroup() %>%
      arrange(desc(gset))
    
    m.all.gset <- ggplot(jfcs.all.annot, aes(x = H3K4me3, y = H3K27me3, color = gset)) + 
      geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle("Hist mod levels of bins, genes split into two groups")
    print(m.all.gset)
    
    m.all.gset.linear <- ggplot(jfcs.all.annot, aes(x = exp(H3K4me3), y = exp(H3K27me3), color = gset)) + 
      geom_point(alpha = 0.2) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle("Hist mod levels of bins, genes split into two groups")
    print(m.all.gset.linear)
    
    # show boxplot
    m.all.gset.dens.act <- ggplot(jfcs.all.annot, aes(x = H3K4me3, fill = gset)) + 
      geom_density(alpha = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle("Hist mod levels of bins, genes split into two groups") + 
      xlab("H3K4me3 levels in HSCs (log scale)")
    print(m.all.gset.dens.act)
    
    m.all.gset.dens.act.lin <- ggplot(jfcs.all.annot, aes(x = exp(H3K4me3), fill = gset)) + 
      geom_density(alpha = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle("Hist mod levels of bins, genes split into two groups") + 
      xlab("H3K4me3 levels in HSCs (linear scale)")
    print(m.all.gset.dens.act.lin)
    
    m.all.gset.dens.repress <- ggplot(jfcs.all.annot, aes(x = H3K27me3, fill = gset)) + 
      geom_density(alpha = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle("Hist mod levels of bins, genes split into two groups") + 
      xlab("H3K27me3 levels in HSCs (log scale)")
    print(m.all.gset.dens.repress)
    
    m.all.gset.dens.repress.lin <- ggplot(jfcs.all.annot, aes(x = exp(H3K27me3), fill = gset)) + 
      geom_density(alpha = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle("Hist mod levels of bins, genes split into two groups") + 
      xlab("H3K27me3 levels in HSCs (linear scale)")
    print(m.all.gset.dens.repress.lin)
    
    
    
  }
  
  
  
 
  
  

  # loop across end points --------------------------------------------------

  gset.spec.str <- paste0(gset.spec, "-specGenes")
  gset.other.str <- paste0(paste(gset.other, collapse = "&"), "-specGenes")

  jbins.spec <- subset(fits.bygenesets.long, geneset %in% gset.spec)$bin
  jbins.other <- subset(fits.bygenesets.long, geneset %in% gset.other)$bin
  
  # remove overlaps
  jbins.overlap <- intersect(jbins.spec, jbins.other)
  
  jbins.spec.filt <- jbins.spec[!jbins.spec %in% jbins.overlap]
  jbins.other.filt <- jbins.other[!jbins.other %in% jbins.overlap]
  
  jints.sub <- subset(jfits.mat.ints, bin %in% c(jbins.spec.filt, jbins.other.filt)) %>%
    rowwise() %>%
    mutate(gset = ifelse(bin %in% jbins.spec.filt, 
                         gset.spec.str, 
                         gset.other.str))
                         # paste0(gset.spec, "-specGenes"), 
                         # paste0(paste(gset.other, collapse = "&"), "-specGenes")))
  jints.sub$gset <- factor(jints.sub$gset, levels = c(gset.spec.str, gset.other.str))
  
  
  jfcs.sub <- subset(jfits.mat.logfcs, bin %in% c(jbins.spec.filt, jbins.other.filt)) %>%
    rowwise() %>%
    mutate(gset = ifelse(bin %in% jbins.spec.filt, 
                         # paste0(gset.spec, "-specGenes"), 
                         # paste0(paste(gset.other, collapse = "&"), "-specGenes"))) %>%
                         gset.spec.str, 
                         gset.other.str)) %>%
    left_join(., subset(jints.sub, select = c(bin, H3K4me1, H3K4me3, H3K27me3)), by = "bin")
  jfcs.sub$gset <- factor(jfcs.sub$gset, levels = c(gset.spec.str, gset.other.str))
  
  
  m.cloud <- ggplot(jints.sub, aes(x = exp(H3K4me3), y = exp(H3K27me3), color = gset)) + 
    geom_point(alpha = 0.1) + 
    facet_wrap(~gset) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle)
  print(m.cloud)
  
  jscale <- 1
  m.cloud.arrows <- ggplot(jfcs.sub, aes(x = exp(H3K4me3), y = exp(H3K27me3), color = gset)) + 
    geom_point(alpha = 0.1) + 
    facet_wrap(~gset) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_segment(mapping = aes(xend = jscale * exp(H3K4me3 + H3K4me3.fc), yend = jscale * exp(H3K27me3 + H3K27me3.fc)),
                 arrow = arrow(length=unit(0.25,"cm"), ends = "last"), alpha = 0.4, size = 0.1) + 
    ggtitle(jtitle)
  print(m.cloud.arrows)
  
  m.cloud.log <- ggplot(jfcs.sub, aes(x = H3K4me3, y = H3K27me3, color = gset)) + 
    geom_point(alpha = 0.1) + 
    facet_wrap(~gset) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle)
  print(m.cloud.log)
  
  jscale2 <- 0.2
  m.cloud.arrows <- ggplot(jfcs.sub, aes(x = H3K4me3, y = H3K27me3, color = gset)) + 
    # geom_point(alpha = 0.1) + 
    facet_wrap(~gset) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_segment(mapping = aes(xend = H3K4me3 + jscale2 * H3K4me3.fc, yend = H3K27me3 + jscale2 * H3K27me3.fc),
                 arrow = arrow(length=unit(0.25,"cm"), ends = "last"), alpha = 0.4, size = 0.1) + 
    ggtitle(jtitle)
  print(m.cloud.arrows)
  
  m.cloud.arrows.zoom <- ggplot(jfcs.sub, aes(x = H3K4me3, y = H3K27me3, color = gset)) + 
    # geom_point(alpha = 0.1) + 
    facet_wrap(~gset) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_segment(mapping = aes(xend = H3K4me3 + jscale2 * H3K4me3.fc, yend = H3K27me3 + jscale2 * H3K27me3.fc),
                 arrow = arrow(length=unit(0.25,"cm"), ends = "last"), alpha = 0.4, size = 0.1) + 
    ggtitle(jtitle) + 
    coord_cartesian(ylim=c(-12, -10), xlim = c(-9, -11)) 
  print(m.cloud.arrows.zoom)
  
  # plot densities to show they are different directions
  
  m.fc.h3k4me1 <- ggplot(jfcs.sub, aes(x = H3K4me1.fc, fill = gset)) + 
    geom_density(alpha = 0.4) + 
    geom_vline(xintercept = 0) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle)
  print(m.fc.h3k4me1)
  
  m.fc.h3k4me3 <- ggplot(jfcs.sub, aes(x = H3K4me3.fc, fill = gset)) + 
    geom_density(alpha = 0.4) + 
    geom_vline(xintercept = 0) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle)
  print(m.fc.h3k4me3)
  
  m.fc.h3k27me3 <- ggplot(jfcs.sub, aes(x = H3K27me3.fc, fill = gset)) + 
    geom_density(alpha = 0.4) + 
    geom_vline(xintercept = 0) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle)
  print(m.fc.h3k27me3)
  
  # show h3k4me3 vs h3k27me3
  nbins1 <- nrow(jfcs.sub)
  # check Bcells-epcific high K27me3 and K4me3 FC
  
  
  if (ctype.end == "ClusterGranulocytes"){
    print("printing exceptions for Granu")
    jcheck <- subset(jfcs.sub, gset == "Neutrophil-specGenes") %>%
      ungroup() %>%
      filter(H3K4me3.fc > 0) %>%
      arrange(H3K27me3.fc) %>%
      mutate(jrank = rank(H3K27me3.fc))
    jbins.lab <- unique(subset(jcheck, jrank <= 25)$bin)
    jfcs.sub.lab <- jfcs.sub %>%
      rowwise() %>%
      mutate(glab = ifelse(bin %in% jbins.lab, gene, NA))
    
    m.fc.k4me3_vs_k27me3.lab <- ggplot(jfcs.sub.lab, aes(x = H3K4me3.fc, y = H3K27me3.fc, color = gset, label = glab)) + 
      geom_point(alpha = 0.4) + 
      geom_text_repel(color = 'black') + 
      geom_vline(xintercept = 0) + 
      geom_hline(yintercept = 0) + 
      facet_wrap(~gset) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle(jtitle, paste("nbins:", nbins1, "print top granu genes"))
    print(m.fc.k4me3_vs_k27me3.lab)
    
  }
  
  
  m.fc.k4me3_vs_k27me3 <- ggplot(jfcs.sub, aes(x = H3K4me3.fc, y = H3K27me3.fc, color = gset)) + 
    geom_point(alpha = 0.4) + 
    geom_density_2d(alpha = 0.8, color = "black") + 
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0) + 
    facet_wrap(~gset) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle, paste("nbins:", nbins1))
  print(m.fc.k4me3_vs_k27me3)
  
  xrange.fc <- ggplot_build(m.fc.k4me3_vs_k27me3)$layout$panel_scales_x[[1]]$range$range
  yrange.fc <- ggplot_build(m.fc.k4me3_vs_k27me3)$layout$panel_scales_y[[1]]$range$range
  
  # add CI
  jfcs.sub.merge <- left_join(jfcs.sub, jmat.fc.ci.lst[[ctype.end]], by = "bin")
    # filter(H3K4me3.fc.upper - H3K27me3.fc.lower <= fc.max,
    #        H3K27me3.fc.upper - H3K27me3.fc.lower <= fc.max)
  nbins2 <- nrow(jfcs.sub.merge)
  
  
  m.fc.k4me3_vs_k27me3.ci <- ggplot(jfcs.sub.merge %>% filter(!is.na(H3K4me3.fc) & !is.na(H3K27me3.fc)), aes(x = H3K4me3.fc, y = H3K27me3.fc, color = gset)) + 
    geom_errorbar(mapping = aes(ymin = H3K27me3.fc.lower, ymax = H3K27me3.fc.upper), alpha = 0.1, width = 0) + 
    geom_errorbarh(mapping = aes(xmin = H3K4me3.fc.lower, xmax = H3K4me3.fc.upper), alpha = 0.1, height = 0) + 
    geom_point(alpha = 0.4) + 
    geom_density_2d(alpha = 0.5, color = "black") + 
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0) + 
    facet_wrap(~gset) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle, paste("Nbins:", nbins2)) + 
    coord_cartesian(xlim = xrange.fc, ylim = yrange.fc)
  print(m.fc.k4me3_vs_k27me3.ci)
  
  if (ctype.end == "ClusterBcells"){
    print("Plotting exceptions for Bcells ")
    # jcheck <- jfcs.sub %>%
    #   arrange(desc(sqrt(H3K4me3.fc ^ 2 + H3K27me3.fc ^ 2)))
    # jfcs.sub.igk <- jfcs.sub %>%
    jfcs.sub.igk <- jfcs.sub.merge %>%
      mutate(is.igkv = ifelse(grepl("Igkv", gene), gene, NA))
    m.fc.k4me3_vs_k27me3.label_igk <-  ggplot(jfcs.sub.igk, aes(x = H3K4me3.fc, y = H3K27me3.fc, color = gset, label = is.igkv)) + 
      geom_point(alpha = 0.4) + 
      ggrepel::geom_text_repel(size = 2) + 
      geom_density_2d(alpha = 0.8, color = "black") + 
      geom_vline(xintercept = 0) + 
      geom_hline(yintercept = 0) + 
      facet_wrap(~gset) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle(jtitle, paste("nbins:", nbins1))
    print(m.fc.k4me3_vs_k27me3.label_igk)
    
    # plot without Igk genes
    m.fc.k4me3_vs_k27me3.no_igk <-  ggplot(jfcs.sub.igk %>% filter(is.na(is.igkv) & !is.na(H3K4me3.fc) & !is.na(H3K27me3.fc)), aes(x = H3K4me3.fc, y = H3K27me3.fc, color = gset)) + 
      geom_point(alpha = 0.4) + 
      geom_density_2d(alpha = 0.8, color = "black") + 
      geom_vline(xintercept = 0) + 
      geom_hline(yintercept = 0) + 
      facet_wrap(~gset) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle(jtitle, paste("nbins:", nbins1, "Igkv genes removed"))
    print(m.fc.k4me3_vs_k27me3.no_igk)
    
    # plot without Igk genes
    m.fc.k4me3_vs_k27me3.no_igk.ci <-  ggplot(jfcs.sub.igk %>% filter(is.na(is.igkv)), aes(x = H3K4me3.fc, y = H3K27me3.fc, color = gset)) + 
      geom_errorbar(mapping = aes(ymin = H3K27me3.fc.lower, ymax = H3K27me3.fc.upper), alpha = 0.1, width = 0) + 
      geom_errorbarh(mapping = aes(xmin = H3K4me3.fc.lower, xmax = H3K4me3.fc.upper), alpha = 0.1, height = 0) + 
      geom_point(alpha = 0.4) + 
      geom_density_2d(alpha = 0.8, color = "black") + 
      geom_vline(xintercept = 0) + 
      geom_hline(yintercept = 0) + 
      facet_wrap(~gset) + 
      theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      ggtitle(jtitle, paste("nbins:", nbins1, "Igkv genes removed"))
    print(m.fc.k4me3_vs_k27me3.no_igk.ci)
    
    
  }
  
  # print(as.data.frame(head(subset(jfcs.sub.merge, gset == "Neutrophil-specGenes" & !is.na(H3K4me3.fc)) %>% arrange(H3K27me3.fc))))
  
})



fewer.k27me3 <- TRUE
jdate <- "2020-06-05"
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.setup_for_poisson_regression"
assertthat::assert_that(dir.exists(indir))
jprefix <- file.path(indir, paste0("integrated_analysis_3_marks.stringentDE.WithHighLowExprs.singlecells.fewerk27me3_", fewer.k27me3, ".forPoissonRegression.CountR1only.", jdate))
infrdata <- paste0(jprefix, ".smaller.RData")
load(infrdata, v=T)


fitsdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits.OtherNorm.DiffDists"
fitprefix <- file.path(fitsdir, paste0("fit_poisson_model_on_TSS.MouseBM.NormMeth_", jnorm, ".bsize_", bsize))
infits <- paste0(fitprefix, ".RData")
load(infits, v=T)


# jgenes <- c("Hlf", "Hbb-bh2", "Pax5", "S100a7a", "Irf4", "Meis1", "Hoxa4")
jgenes <- c("Hlf", "Tead1", "Hoxa9", "Meis1", "Hoxb9", "Hoxb3", "Hoxd4", "Adgrg1", "S100a7a", "S100a8", "Ltf", "Hdac4", "Hbb-bs", "Hbb-y", "Sox6", "Ebf1", "Irf4", "Pax5", "Cd180", "Cd38", "Iglv3", "Bach2")
for (jgene in jgenes){
  for (jmark in jmarks){
    (jbin <- rownames(tss.mats.filt.fromref.cellfilt[[jmark]])[grepl(jgene, rownames(tss.mats.filt.fromref.cellfilt[[jmark]]))])
    jrow <- tss.mats.filt.fromref.cellfilt[[jmark]][jbin, ]
    cnames <- colnames(tss.mats.filt.fromref.cellfilt[[jmark]])
    dat.annots.filt.mark <- dat.annots.filt.forfit[[jmark]]
    ncuts.cells.mark <- ncuts.cells[[jmark]]
    refit <- RefitPoissonForPlot(jrow = jrow, cnames = cnames, dat.annots.filt.mark = dat.annots.filt.mark, ncuts.cells.mark = ncuts.cells.mark)
    
    m.raw.log <- ggplot(refit$input.dat, aes(x = Cluster, y = logLambda)) + 
      geom_errorbar(mapping = aes(ymin = logLambdaLower, ymax = logLambdaUpper), data = refit$params.mean.dat, width = 0.1, color = 'white') + 
      geom_jitter(width = 0.1, height = 0, alpha = 0.1) + 
      # geom_hline(yintercept = refit$params.int.dat$logLambda, linetype = "dotted") + 
      theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      xlab("" ) + ylab("log(ncuts) - log(total cuts)") + 
      ggtitle(paste(jgene, jbin, jmark)) 
    print(m.raw.log)
    
    
    m.fit <- ggplot(refit$input.dat, aes(x = Cluster, y = logLambda)) + 
      geom_jitter(width = 0.1, height = 0, alpha = 0.1) + 
      geom_errorbar(mapping = aes(ymin = logLambdaLower, ymax = logLambdaUpper), data = refit$params.mean.dat, width = 0.1) + 
      geom_hline(yintercept = refit$params.int.dat$logLambda, linetype = "dotted") + 
      theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      xlab("" ) + ylab("log(ncuts) - log(total cuts)") + 
      ggtitle(paste(jgene, jbin, jmark)) 
    print(m.fit)
    
    # fit in linear scale
    
    m.raw.linear <- ggplot(refit$input.dat, aes(x = Cluster, y = exp(logLambda))) + 
      geom_jitter(width = 0.1, height = 0, alpha = 0.1) + 
      # geom_hline(yintercept = exp(refit$params.int.dat$logLambda), linetype = "dotted") + 
      theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      xlab("" ) + ylab("ncuts / totalcuts") + 
      ggtitle(paste(jgene, jbin, jmark)) 
    print(m.raw.linear)
    
    m.fit.linear <- ggplot(refit$input.dat, aes(x = Cluster, y = exp(logLambda))) + 
      geom_jitter(width = 0.1, height = 0, alpha = 0.1) + 
      geom_errorbar(mapping = aes(ymin = exp(logLambdaLower), ymax = exp(logLambdaUpper)), 
                    data = refit$params.mean.dat, width = 0.1) + 
      geom_hline(yintercept = exp(refit$params.int.dat$logLambda), linetype = "dotted") + 
      theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      xlab("" ) + ylab("ncuts / totalcuts") + 
      ggtitle(paste(jgene, jbin, jmark)) 
    print(m.fit.linear)
  }
}


if (make.plots){
  dev.off()
}






 
