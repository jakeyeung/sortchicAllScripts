# Jake Yeung
# Date of Creation: 2021-02-03
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_deeper/2-compare_logFC_bins_HSPC_vs_nonHSPCs.R
# 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

jstart <- Sys.time()

make.plots <- TRUE


# pvalcutoff <- 1e-100
padjcutoff <- 1e-50
# jfactor <- 1.96  # 95% confidence interval
jfactor <- 2.576  # 99% confidence interval
jsize <- 0.05

# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
  

indir.fits <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.HSPCs_vs_nonHSPCs"

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/compare_DE_analysis_spikein_vs_total"
pdfname <- paste0("HSPCs_vs_nonHSPCs.total_vs_spikeins_plots.with_se.", jfactor, ".", Sys.Date(), ".pdf")
outpdf <- file.path(outdir, pdfname)


# Load DEs spikeins vs normals  -------------------------------------------

jsuffix.spikeins <- "HSPCs_vs_nonHSPCs.spikeins"
jdate.spikeins <- "2020-12-29"
# jsuffix.total <- "total"
jsuffix.total <- "HSPCs_vs_nonHSPCs.total"
jdate.total <- "2021-02-03"
out.lst.spikeins <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0("poisson_fit_bins.", jmark, ".", jdate.spikeins, ".newannot2.witherrors.MoreBins.", jsuffix.spikeins, ".RData")
  inf.fits <- file.path(indir.fits, fname)
  print(inf.fits)
  load(inf.fits, v=T)
  params.long <- SummarizeParamsPvalues(jfits.lst, jmark = jmark, paramname = "Cluster") %>%
    mutate(log2fc = estimate / log(2),
           estimate = estimate / log(2), 
           se = jfactor * se / log(2))
  params.long$padj <- p.adjust(params.long$pval.param)
  jnames <- names(jfits.lst); names(jnames) <- jnames
  pvals.long <- lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xvec <- x$pval
    data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  return(list(params.long = params.long, pvals.long = pvals.long))
})


out.lst.total <- lapply(jmarks, function(jmark){
  print(jmark)
  fname <- paste0("poisson_fit_bins.", jmark, ".", jdate.total, ".newannot2.witherrors.MoreBins.", jsuffix.total, ".RData")
  inf.fits <- file.path(indir.fits, fname)
  load(inf.fits, v=T)
  params.long <- SummarizeParamsPvalues(jfits.lst, jmark = jmark, paramname = "Cluster") %>%
    mutate(log2fc = estimate / log(2),
           estimate = estimate / log(2), 
           se = jfactor * se / log(2))
  params.long$padj <- p.adjust(params.long$pval.param)
  jnames <- names(jfits.lst); names(jnames) <- jnames
  pvals.long <- lapply(jnames, function(jname){
    x <- jfits.lst[[jname]]
    xvec <- x$pval
    data.frame(bin = jname, pval = xvec, mark = jmark, stringsAsFactors = FALSE)
  }) %>%
    bind_rows()
  return(list(params.long = params.long, pvals.long = pvals.long))
})


print(head(out.lst.total$H3K4me1$params.long))
print(head(out.lst.spikeins$H3K4me1$params.long))


params.spikeins.long <- lapply(out.lst.spikeins, function(out){
  out$params.long
}) %>%
  bind_rows()

pvals.spikeins.long <- lapply(out.lst.spikeins, function(out){
  out$pvals.long
}) %>%
  bind_rows()


params.total.long <- lapply(out.lst.total, function(out){
  out$params.long
}) %>%
  bind_rows()

pvals.total.long <- lapply(out.lst.total, function(out){
  out$pvals.long
}) %>%
  bind_rows()



# Check distributions  ----------------------------------------------------

if (make.plots){
  pdf(outpdf, useDingbats = FALSE)
}

ggplot(params.total.long, aes(x = log2fc, fill = mark)) + 
  geom_density(alpha = 0.25) + 
  facet_wrap(~mark, nrow = 1) + 
  coord_cartesian(c(-5, 5)) + 
  ggtitle("Norm total cuts") + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(params.spikeins.long, aes(x = log2fc, fill = mark)) + 
  geom_density(alpha = 0.25) + 
  facet_wrap(~mark, nrow = 1) + 
  coord_cartesian(c(-5, 5)) + 
  ggtitle("Norm spikein cuts") + 
  geom_vline(xintercept = 0, linetype = "dotted") + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")


# Merge and plot log2FC  ---------------------------------------------------

params.merge.lst <- lapply(jmarks, function(jmark){
  jsub.total <- subset(params.total.long, mark == jmark) %>%
    dplyr::rename(estimate.total = estimate, se.total = se, z.total = z, pval.param.total = pval.param, log2fc.total = log2fc, padj.total = padj)
  jsub.spikeins <- subset(params.spikeins.long, mark == jmark) %>%
    dplyr::rename(estimate.spikeins = estimate, se.spikeins = se, z.spikeins = z, pval.param.spikeins = pval.param, log2fc.spikeins = log2fc, padj.spikeins = padj)
  jsub.merge <- left_join(jsub.total, jsub.spikeins, by = c("bin", "param", "mark"))
  return(jsub.merge)
})

m.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(params.merge.lst[[jmark]] %>% filter(abs(estimate.total) < 5 & abs(estimate.spikeins) < 5), aes(x = estimate.total, y = estimate.spikeins)) + 
    geom_point(alpha = 0.25) + 
    # geom_errorbarh(mapping = aes(xmin = estimate.total - se.total, xmax = estimate.total + se.total), alpha = 0.05) + 
    # geom_errorbar(mapping = aes(ymin = estimate.spikeins - se.spikeins, ymax = estimate.spikeins + se.spikeins), alpha = 0.05) + 
    geom_density_2d() + 
    theme_bw() + 
    ggtitle(jmark) + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    # coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

JFuncs::multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 2)

# add error bars
m.lst <- lapply(jmarks, function(jmark){
  
  jsub <- params.merge.lst[[jmark]] %>% filter(abs(estimate.total) < 5 & abs(estimate.spikeins) < 5) %>%
    filter(se.total < 10)
  
  m1 <- ggplot(jsub, aes(x = estimate.total, y = se.total)) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + 
    ggtitle(paste(jmark, "all points")) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m1)
  
  m1 <- ggplot(jsub, aes(x = se.total)) + 
    geom_histogram(bins = 100) + 
    theme_bw() + 
    ggtitle(paste(jmark, "all points")) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m1)
  
  # plot(jsub$estimate.total, jsub$se.total, pch = 20)
  # plot(jsub$estimate.spikeins, jsub$se.spikeins, pch = 20)
  # plot(density(jsub$se.total))
  
  m <- ggplot(jsub %>% sample_frac(., size = jsize), aes(x = estimate.total, y = estimate.spikeins)) + 
    geom_point(alpha = 0.25) + 
    geom_errorbarh(mapping = aes(xmin = estimate.total - se.total, xmax = estimate.total + se.total), alpha = 0.05) + 
    geom_errorbar(mapping = aes(ymin = estimate.spikeins - se.spikeins, ymax = estimate.spikeins + se.spikeins), alpha = 0.05) + 
    # geom_density_2d() + 
    theme_bw() + 
    ggtitle(paste(jmark, "subsampled", jsize)) + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

JFuncs::multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 2)



# Get bins where signs are flipped  ---------------------------------------

bins.flipped.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  jfilt <- subset(params.merge.lst[[jmark]], estimate.total - se.total > 0 & estimate.spikeins + se.spikeins < 0)
  return(jfilt$bin)
})



# Compare between marks  --------------------------------------------------

print(head(params.total.long))

params.total.wide <- data.table::dcast(data = params.total.long, formula = bin ~ mark, value.var = "estimate", fill = NA)
params.total.wide <- params.total.wide[complete.cases(params.total.wide), ]

params.total.wide.se <- data.table::dcast(data = params.total.long, formula = bin ~ mark, value.var = "se", fill = NA)
params.total.wide.se <- params.total.wide.se[complete.cases(params.total.wide.se), ]
assertthat::assert_that(identical(params.total.wide$bin, params.total.wide.se$bin))

params.total.wide.lower <- as.data.frame(as.matrix(params.total.wide[, -1]) - as.matrix(params.total.wide.se[, -1]))
params.total.wide.upper <- as.data.frame(as.matrix(params.total.wide[, -1]) + as.matrix(params.total.wide.se[, -1]))
rownames(params.total.wide.lower) <- params.total.wide$bin
rownames(params.total.wide.upper) <- params.total.wide$bin

colnames(params.total.wide.lower) <- paste(colnames(params.total.wide.lower), "lower", sep = ".")
colnames(params.total.wide.upper) <- paste(colnames(params.total.wide.upper), "upper", sep = ".")

# merge with estimates with se into a wide DF
params.total.wide <- left_join(params.total.wide, 
                               data.frame(bin = rownames(params.total.wide.lower), params.total.wide.lower, stringsAsFactors = FALSE))
params.total.wide <- left_join(params.total.wide, 
                               data.frame(bin = rownames(params.total.wide.upper), params.total.wide.upper, stringsAsFactors = FALSE))



  
params.spikeins.wide <- data.table::dcast(data = params.spikeins.long, formula = bin ~ mark, value.var = "estimate", fill = NA)
params.spikeins.wide <- params.spikeins.wide[complete.cases(params.spikeins.wide), ]

params.spikeins.wide.se <- data.table::dcast(data = params.spikeins.long, formula = bin ~ mark, value.var = "se", fill = NA)
params.spikeins.wide.se <- params.spikeins.wide.se[complete.cases(params.spikeins.wide.se), ]
assertthat::assert_that(identical(params.spikeins.wide$bin, params.spikeins.wide.se$bin))

params.spikeins.wide.lower <- as.data.frame(as.matrix(params.spikeins.wide[, -1]) - as.matrix(params.spikeins.wide.se[, -1]))
params.spikeins.wide.upper <- as.data.frame(as.matrix(params.spikeins.wide[, -1]) + as.matrix(params.spikeins.wide.se[, -1]))
rownames(params.spikeins.wide.lower) <- params.spikeins.wide$bin
rownames(params.spikeins.wide.upper) <- params.spikeins.wide$bin

colnames(params.spikeins.wide.lower) <- paste(colnames(params.spikeins.wide.lower), "lower", sep = ".")
colnames(params.spikeins.wide.upper) <- paste(colnames(params.spikeins.wide.upper), "upper", sep = ".")



# merge with estimates with se into a wide DF
params.spikeins.wide <- left_join(params.spikeins.wide, 
                               data.frame(bin = rownames(params.spikeins.wide.lower), params.spikeins.wide.lower, stringsAsFactors = FALSE))
params.spikeins.wide <- left_join(params.spikeins.wide, 
                               data.frame(bin = rownames(params.spikeins.wide.upper), params.spikeins.wide.upper, stringsAsFactors = FALSE))

jmarkref <- "H3K27me3"
for (jmark in jmarks[jmarks != jmarkref]){
  m <- ggplot(params.total.wide, aes_string(x = jmarkref, y = jmark)) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + 
    geom_density_2d() + 
    ggtitle("log2FC norm by total") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(params.spikeins.wide, aes_string(x = jmarkref, y = jmark)) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + 
    geom_density_2d() + 
    ggtitle("log2FC norm by spikeins") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
}

# get signif bins 
jfilt <- subset(params.spikeins.long, mark == jmarkref) %>%
  arrange(estimate) %>%
  filter(abs(estimate) < 5 & pval.param < padjcutoff & estimate < 0)

bins.lost <- unique(jfilt$bin)
print(length(bins.lost))



# Replot lost bins --------------------------------------------------------

jmarkref <- "H3K27me3"
for (jmark in jmarks[jmarks != jmarkref]){
  m <- ggplot(params.total.wide %>% filter(bin %in% bins.lost), aes_string(x = jmarkref, y = jmark)) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + 
    geom_density_2d() + 
    ggtitle("log2FC norm by total bins lost only") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m2 <- m + geom_errorbarh(mapping = aes_string(xmin = paste0(jmarkref, ".lower"), xmax = paste0(jmarkref, ".upper")), alpha = 0.25) + 
    geom_errorbar(mapping = aes_string(ymin = paste0(jmark, ".lower"), ymax = paste0(jmark, ".upper")), alpha = 0.25) 
  print(m2)
  
  m <- ggplot(params.spikeins.wide %>% filter(bin %in% bins.lost), aes_string(x = jmarkref, y = jmark)) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + 
    geom_density_2d() + 
    ggtitle("log2FC norm by spikeins bins lost only") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m2 <- m + geom_errorbarh(mapping = aes_string(xmin = paste0(jmarkref, ".lower"), xmax = paste0(jmarkref, ".upper")), alpha = 0.25) + 
    geom_errorbar(mapping = aes_string(ymin = paste0(jmark, ".lower"), ymax = paste0(jmark, ".upper")), alpha = 0.25) 
  print(m2)
  
}

# JFuncs::multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 2)
m.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  m <- ggplot(params.merge.lst[[jmark]] %>% filter(bin %in% bins.lost & abs(estimate.total) < 5 & abs(estimate.spikeins) < 5), aes(x = estimate.total, y = estimate.spikeins)) +
    geom_point(alpha = 0.25) + 
    geom_density_2d() + 
    theme_bw() + 
    ggtitle(jmark, "bins lost only") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m2 <- m + geom_errorbarh(mapping = aes(xmin = estimate.total - se.total, xmax = estimate.total + se.total), alpha = 0.25) + 
    geom_errorbar(mapping = aes(ymin = estimate.spikeins - se.spikeins, ymax = estimate.spikeins + se.spikeins), alpha = 0.25) 
  print(m2)
  return(m)
})
JFuncs::multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 2)


# Replot flipped bins  ----------------------------------------------------



jmarkref <- "H3K27me3"
for (jmark in jmarks[jmarks != jmarkref]){
  m <- ggplot(params.total.wide %>% filter(bin %in% bins.flipped.lst[[jmarkref]]), aes_string(x = jmarkref, y = jmark)) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + 
    geom_density_2d() + 
    ggtitle("log2FC norm by total bins flipped only") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  m2 <- m + geom_errorbarh(mapping = aes_string(xmin = paste0(jmarkref, ".lower"), xmax = paste0(jmarkref, ".upper")), alpha = 0.25) + 
    geom_errorbar(mapping = aes_string(ymin = paste0(jmark, ".lower"), ymax = paste0(jmark, ".upper")), alpha = 0.25) 
  print(m2)
  
  m <- ggplot(params.spikeins.wide %>% filter(bin %in% bins.flipped.lst[[jmarkref]]), aes_string(x = jmarkref, y = jmark)) + 
    geom_point(alpha = 0.25) + 
    theme_bw() + 
    geom_density_2d() + 
    ggtitle("log2FC norm by spikeins bins flipped only") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  m2 <- m + geom_errorbarh(mapping = aes_string(xmin = paste0(jmarkref, ".lower"), xmax = paste0(jmarkref, ".upper")), alpha = 0.25) + 
    geom_errorbar(mapping = aes_string(ymin = paste0(jmark, ".lower"), ymax = paste0(jmark, ".upper")), alpha = 0.25) 
  print(m2)
  print(m)
}

# JFuncs::multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 2)
m.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  m <- ggplot(params.merge.lst[[jmark]] %>% filter(bin %in% bins.flipped.lst[[jmark]] & abs(estimate.total) < 5 & abs(estimate.spikeins) < 5), aes(x = estimate.total, y = estimate.spikeins)) +
    geom_point(alpha = 0.25) + 
    geom_density_2d() + 
    theme_bw() + 
    ggtitle(jmark, "bins flipped only") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
  m2 <- m + geom_errorbarh(mapping = aes(xmin = estimate.total - se.total, xmax = estimate.total + se.total), alpha = 0.25) + 
    geom_errorbar(mapping = aes(ymin = estimate.spikeins - se.spikeins, ymax = estimate.spikeins + se.spikeins), alpha = 0.25) 
  
  print(m)
  print(m2)
  
  return(m)
})
JFuncs::multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 2)


# Write filpped bins but colored ------------------------------------------


m.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  
  jsub <- params.merge.lst[[jmark]] %>% 
    filter(abs(estimate.total) < 5 & abs(estimate.spikeins) < 5) %>% 
    mutate(is.flipped = bin %in% bins.flipped.lst[[jmark]]) %>% 
    arrange(is.flipped)
  
  m <- ggplot(jsub, 
              aes(x = estimate.total, y = estimate.spikeins, color = is.flipped)) +
    geom_point(alpha = 1) + 
    theme_bw() + 
    ggtitle(jmark, "bins flipped only") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
  
  m2 <- m + geom_errorbarh(mapping = aes(xmin = estimate.total - se.total, xmax = estimate.total + se.total), alpha = 0.25) + 
    geom_errorbar(mapping = aes(ymin = estimate.spikeins - se.spikeins, ymax = estimate.spikeins + se.spikeins), alpha = 0.25) 
  
  print(m)
  print(m2)
  
  return(m)
})
JFuncs::multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 2)


# flip axes
m.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  
  jsub <- params.merge.lst[[jmark]] %>% 
    filter(abs(estimate.total) < 5 & abs(estimate.spikeins) < 5) %>% 
    mutate(is.flipped = bin %in% bins.flipped.lst[[jmark]]) %>% 
    arrange(is.flipped)
  
  m <- ggplot(jsub, 
              aes(y = estimate.total, x = estimate.spikeins, color = is.flipped)) +
    geom_point(alpha = 1) + 
    theme_bw() + 
    ggtitle(jmark, "bins flipped only") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
  
  m2 <- m + geom_errorbar(mapping = aes(ymin = estimate.total - se.total, ymax = estimate.total + se.total), alpha = 0.25) + 
    geom_errorbarh(mapping = aes(xmin = estimate.spikeins - se.spikeins, xmax = estimate.spikeins + se.spikeins), alpha = 0.25) 
  
  print(m)
  print(m2)
  
  return(m)
})
JFuncs::multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 2)





if (make.plots){
  dev.off()
}


# Write flipped bins to output -------------------------------------------

dat.bins.flipped.lst <- lapply(bins.flipped.lst, function(x){
  GetBedFromCoords(x, add.chr = FALSE, strip.chr = FALSE)
})

dat.bins.lost <- GetBedFromCoords(bins.lost)

for (jmark in jmarks){
  print(jmark)
  outname <- paste0("bins_flipped.", jmark, ".", jfactor, ".", Sys.Date(), ".txt")
  outftmp <- file.path(outdir, outname)
  fwrite(x = dat.bins.flipped.lst[[jmark]], file = outftmp, sep = "\t")
}

outnamelost <- paste0("bins_lost.H3K27me3.", jfactor, ".", Sys.Date(), ".txt")
outftmplost <- file.path(outdir, outnamelost)
fwrite(dat.bins.lost, file = outftmplost, sep = "\t")



# Write objects -----------------------------------------------------------

outfmerged <- file.path(outdir, paste0("params_merge_lst.", Sys.Date(), ".rds"))
saveRDS(params.merge.lst, file = outfmerged)

outftotal <- file.path(outdir, paste0("params_total_wide.", Sys.Date(), ".rds"))
saveRDS(params.total.wide, file = outftotal)

outfspikeins <- file.path(outdir, paste0("params_spikeins_wide.", Sys.Date(), ".rds"))
saveRDS(params.spikeins.wide, file = outfspikeins)





print(Sys.time() - jstart)



# Show nonflipped bins but quadrant 2  ----------------------------------------------------


m.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  
  jsub <- params.merge.lst[[jmark]] %>% 
    filter(abs(estimate.total) < 5 & abs(estimate.spikeins) < 5) %>% 
    mutate(is.flipped = bin %in% bins.flipped.lst[[jmark]]) %>% 
    arrange(is.flipped) %>%
    filter(estimate.spikeins < 0 & estimate.total > 0)
  
  m <- ggplot(jsub %>% filter(!is.flipped), 
              aes(y = estimate.total, x = estimate.spikeins, color = is.flipped)) +
    geom_point(alpha = 1) + 
    theme_bw() + 
    ggtitle(jmark, "bins flipped only") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
  
  m2 <- m + geom_errorbar(mapping = aes(ymin = estimate.total - se.total, ymax = estimate.total + se.total), alpha = 0.25) + 
    geom_errorbarh(mapping = aes(xmin = estimate.spikeins - se.spikeins, xmax = estimate.spikeins + se.spikeins), alpha = 0.25) 
  
  print(m)
  print(m2)
  
  mcheck <- ggplot(jsub %>% filter(!is.flipped), 
                   aes(y = estimate.total - se.total, x = estimate.spikeins + se.spikeins, color = is.flipped)) +
    geom_point(alpha = 1) + 
    theme_bw() + 
    ggtitle(jmark, "bins flipped only") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  
  return(m)
})
JFuncs::multiplot(m.lst[[1]], m.lst[[2]], m.lst[[3]], m.lst[[4]], cols = 2)





