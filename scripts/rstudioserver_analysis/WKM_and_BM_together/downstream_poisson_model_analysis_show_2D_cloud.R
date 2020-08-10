# Jake Yeung
# Date of Creation: 2020-06-06
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/downstream_poisson_model_analysis_show_2D_cloud.R
# Show 2D cloud

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)

jdate <- "2020-06-06"
jdate2 <- "2020-06-08"

fewer.k27me3 <- TRUE
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.setup_for_poisson_regression"
jprefix <- file.path(indir, paste0("integrated_analysis_3_marks.stringentDE.WithHighLowExprs.singlecells.fewerk27me3_", fewer.k27me3, ".forPoissonRegression.CountR1only.", jdate))

# outfits <- file.path(indir, paste0("fit_poisson_model_on_TSS.", jdate, ".RData"))
infit.wrangled <- file.path(indir, paste0("fit_poisson_model_on_TSS.DownstreamWrangled.", jdate, ".RData"))
pdfout <- file.path(indir, paste0("fit_poisson_model_on_TSS.Downstream2DClouds.", Sys.Date(), ".pdf"))

load(infit.wrangled, v=T)

jfits.long <- jfits.long %>%
  group_by(bin) %>%
  filter(abs(logLambda) < 10) %>%
  filter(abs(logintercept) < 20)

# add arrows? check for Bcell-specific genes vs Eryth+Granu specific genes
print(unique(fits.bygenesets.long$geneset))
print(unique(fits.bygenesets.long$cluster))

# Make into matrix and plot 2D matrix  ------------------------------------


# Plot the 2D cloud -------------------------------------------------------


ctypes.end <- list("ClusterGranulocytes", "ClusterBcells", "ClusterErythroblasts"); names(ctypes.end) <- ctypes.end
gset.specs <- list("Neutrophil", "Bcell", "Erythroblast"); names(gset.specs) <- ctypes.end
gset.others <- list(c("Bcell", "Erythroblast"), c("Erythroblast", "Neutrophil"), c("Bcell", "Neutrophil")); names(gset.others) <- ctypes.end


ctype.end <- ctypes.end[[1]] 

pdf(pdfout, width = 1020/72, height = 815/72, useDingbats = FALSE)

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
  
  if (ctype.end == ctypes.end[[1]]){
    
    
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
    geom_segment(mapping = aes(xend = jscale * exp(H3K4me3 + H3K4me1.fc), yend = jscale * exp(H3K27me3 + H3K27me3.fc)),
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
    geom_segment(mapping = aes(xend = H3K4me3 + jscale2 * H3K4me1.fc, yend = H3K27me3 + jscale2 * H3K27me3.fc),
                 arrow = arrow(length=unit(0.25,"cm"), ends = "last"), alpha = 0.4, size = 0.1) + 
    ggtitle(jtitle)
  print(m.cloud.arrows)
  
  m.cloud.arrows.zoom <- ggplot(jfcs.sub, aes(x = H3K4me3, y = H3K27me3, color = gset)) + 
    # geom_point(alpha = 0.1) + 
    facet_wrap(~gset) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_segment(mapping = aes(xend = H3K4me3 + jscale2 * H3K4me1.fc, yend = H3K27me3 + jscale2 * H3K27me3.fc),
                 arrow = arrow(length=unit(0.25,"cm"), ends = "last"), alpha = 0.4, size = 0.1) + 
    ggtitle(jtitle) + 
    coord_cartesian(ylim=c(-12, -10), xlim = c(-9, -11)) 
  print(m.cloud.arrows.zoom)
  
  # plot densities to show they are different directions
  
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
  m.fc.k4me3_vs_k27me3 <- ggplot(jfcs.sub, aes(x = H3K4me3.fc, y = H3K27me3.fc, color = gset)) + 
    geom_point(alpha = 0.4) + 
    geom_density_2d(alpha = 0.8, color = "black") + 
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0) + 
    facet_wrap(~gset) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jtitle)
  print(m.fc.k4me3_vs_k27me3)
  
})

dev.off()






 
