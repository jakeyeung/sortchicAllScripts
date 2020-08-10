# Jake Yeung
# Date of Creation: 2020-06-05
# File: ~/projects/scchic/scripts/rstudioserver_analysis/WKM_and_BM_together/downstream_poisson_model_analysis.R
# Downstream analysis of poisson model

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(ggrastr)
library(scchicFuncs)
# library(ggrastr)
# library(ggrastr)

# library(devtools)
# dev_mode(TRUE)
# devtools::install_local(path = "/home/jyeung/projects/fromgithub/ggrastr_withjitter")
# library(ggrastr)
# geom_jitter = getFromNamespace("geom_jitter", "ggrastr")

# Functions ---------------------------------------------------------------



# Constants ---------------------------------------------------------------


# jnorm <- "ncuts.alltss"
jnorm <- "ncuts.inbins"
bsize <- 10000

jdate <- "2020-06-05"

fewer.k27me3 <- TRUE
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.setup_for_poisson_regression"
assertthat::assert_that(dir.exists(indir))

jprefix <- file.path(indir, paste0("integrated_analysis_3_marks.stringentDE.WithHighLowExprs.singlecells.fewerk27me3_", fewer.k27me3, ".forPoissonRegression.CountR1only.", jdate))


fitsdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/poisson_fits.OtherNorm.DiffDists"

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/WKM_BM_merged/from_rstudioserver/integrated_analysis_poisson_and_2D_clouds"
assertthat::assert_that(dir.exists(outdir))

fitprefix <- file.path(fitsdir, paste0("fit_poisson_model_on_TSS.MouseBM.NormMeth_", jnorm, ".bsize_", bsize))

infits <- paste0(fitprefix, ".RData")
outfits.wrangled <- paste0(fitprefix, ".DownstreamWrangled.RData")
outpdf <- paste0(fitprefix, ".DownstreamWrangled.pdf")

assertthat::assert_that(!file.exists(outfits.wrangled))
assertthat::assert_that(!file.exists(outpdf))


infrdata <- paste0(jprefix, ".smaller.RData")

assertthat::assert_that(file.exists(infrdata))

load(infrdata, v=T)
load(infits, v=T)

pdf(outpdf, width = 1020/72, height = 815/72, useDingbats = FALSE)

# Wrangle -----------------------------------------------------------------



jfits.dat.lst <- lapply(jmarks, function(jmark){
  jfit.dat <- do.call(rbind, jfits.lst.bymark[[jmark]])
  jfit.dat$bin <- rownames(jfit.dat)
  jfit.dat$gene <- sapply(as.character(jfit.dat$bin), function(b) strsplit(b, ";")[[1]][[2]])
  jfit.dat$ens <- sapply(jfit.dat$gene, function(g) AssignHash(g, g2e, null.fill = g))
  return(jfit.dat)
})

jfits.long.lst <- lapply(jmarks, function(jmark){
  jlong <- jfits.dat.lst[[jmark]] %>%
    dplyr::select(-c(dev.diff, df.diff)) %>%
    reshape2::melt(., id.vars = c("bin", "gene", "ens", "pval"), variable.name = "cluster", value.name = "logLambda")
    # mutate(cluster = ifelse(cluster == "X.Intercept.", "ClusterHSPCs", as.character(cluster)))
  jlong$mark <- jmark
  # extract the X intercept as separate column 
  jlong.noint <- subset(jlong, cluster != "X.Intercept.")
  jlong.int <- subset(jlong, cluster == "X.Intercept.")  %>%
    dplyr::rename(logintercept = logLambda) %>%
    dplyr::select(bin, logintercept, mark)
  jlong.merge <- left_join(jlong.noint, jlong.int)
  return(jlong.merge)
})


# Summarize genoemwide ----------------------------------------------------

head(jfits.long.lst$H3K27me3)

jmark <- "H3K4me3"

# jgene <- "Ebf1"
# (jsub <- subset(jfits.long.lst[[jmark]], gene == jgene))
# ggplot(jsub %>% filter(cluster != "X.Intercept."), aes(x = cluster, y = logLambda)) + geom_point_rast() + 
#   theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  +
#   geom_hline(yintercept = 0) + 
#   ggtitle(paste(jgene, jmark))

m.summary.lst <- lapply(jmarks, function(jmark){
  m.summary <- ggplot(jfits.long.lst[[jmark]], aes(x = cluster, y = logLambda, color = logintercept)) + geom_violin() + 
    # geom_jitter(width = 0.1, height = 0, alpha = 0.1) + 
    scale_color_viridis_c() + 
    theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    ggtitle(jmark, paste("Nbins:", length(unique(jfits.long.lst[[jmark]]$bin))))
})
print(m.summary.lst)

m.summary.boxplot.lst <- lapply(jmarks, function(jmark){
  m.summary <- ggplot(jfits.long.lst[[jmark]], aes(x = cluster, y = logLambda, color = logintercept)) + 
    geom_boxplot() + 
    # geom_jitter(width = 0.1, height = 0, alpha = 0.1) + 
    scale_color_viridis_c() + 
    theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    ggtitle(jmark, paste("Nbins:", length(unique(jfits.long.lst[[jmark]]$bin))))
})
print(m.summary.boxplot.lst)

m.summary.filt.lst <- lapply(jmarks, function(jmark){
  m.summary <- ggplot(jfits.long.lst[[jmark]] %>% filter(abs(logLambda) < 10), aes(x = cluster, y = logLambda, color = logintercept)) + geom_violin() + 
    # geom_jitter(width = 0.1, height = 0, alpha = 0.1) + 
    scale_color_viridis_c() + 
    theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    ggtitle(jmark, paste("Nbins:", length(unique(jfits.long.lst[[jmark]]$bin))))
})
print(m.summary.filt.lst)


m.summary.filt.boxplot.lst <- lapply(jmarks, function(jmark){
  m.summary <- ggplot(jfits.long.lst[[jmark]] %>% filter(abs(logLambda) < 10), aes(x = cluster, y = logLambda, color = logintercept)) + 
    geom_boxplot() + 
    # geom_jitter(width = 0.1, height = 0, alpha = 0.1) + 
    scale_color_viridis_c() + 
    theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    ggtitle(jmark, paste("Nbins:", length(unique(jfits.long.lst[[jmark]]$bin))))
})
print(m.summary.filt.boxplot.lst)



m.volcano.lst <- lapply(jmarks, function(jmark){
  m.summary <- ggplot(jfits.long.lst[[jmark]] %>% filter(abs(logLambda) < 10), aes(x = logLambda, y = -log10(pval), color = logintercept)) +
    geom_point_rast(alpha = 0.05, width = 0, height = 0, size = 6) + 
    scale_color_viridis_c() + 
    facet_wrap(~cluster) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
    theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(jmark, paste("Nbins:", length(unique(jfits.long.lst[[jmark]]$bin))))
})
print(m.volcano.lst)

print(m.volcano.lst$H3K4me3)
print(m.volcano.lst$H3K27me3)
print(m.volcano.lst$H3K4me1)


# # check outliers?
# jsub.high <- subset(jfits.long.lst$H3K4me3, abs(logLambda) < 10) %>%
#   arrange(desc(abs(logLambda)))
# 
# jbin <- sample(jsub.high$bin, size = 1)
# subset(jfits.long.lst$H3K4me3, bin == jbin)
# 
# # plot raw
# counts.dat <- data.frame(cell = colnames(tss.mats.filt.fromref.cellfilt$H3K4me3), ncuts = tss.mats.filt.fromref.cellfilt$H3K4me3[jbin, ], stringsAsFactors = FALSE) %>%
#   left_join(., dat.annots.filt.forfit$H3K4me3) %>%
#   left_join(., ncuts.cells$H3K4me3)
# 
# ggplot(counts.dat, aes(x = Cluster, y = ncuts/ncuts.total)) + geom_jitter(width = 0.1, height = 0) + 
#   theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   ggtitle(jbin)

# Make into long ----------------------------------------------------------

# remove the intercept to show logFoldChanges only 
jfits.long <- lapply(jmarks, function(jmark){
  jtmp <- jfits.long.lst[[jmark]] 
  return(jtmp)
}) %>%
  bind_rows()

jfits.long$mark <- factor(jfits.long$mark, levels = c("H3K4me1", "H3K4me3", "H3K27me3"))


# Summarize by gene sets --------------------------------------------------

genesets <- names(de.ens.sorted.stringent)
names(genesets) <- genesets

fits.bygenesets.long <- lapply(genesets, function(gset){
  ens.keep <- as.character(de.ens.sorted.stringent[[gset]])
  jfits.sub <- subset(jfits.long, ens %in% ens.keep)
  jfits.sub$geneset <- gset
  return(jfits.sub)
}) %>%
  bind_rows() 

ggplot(fits.bygenesets.long %>% filter(abs(logLambda) < 10), aes(x = mark, y = logLambda, fill = cluster)) + 
  geom_boxplot() + theme_bw(24) +
  facet_wrap(~geneset) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  ggtitle("estimated logFCs across gene sets")

m.vol.gsets <- lapply(jmarks, function(jmark){
  m <- ggplot(fits.bygenesets.long %>% filter(abs(logLambda) < 10 & mark == jmark), aes(x = logLambda, y = -log10(pval), color = logintercept)) +
    geom_point_rast(alpha = 0.5, size = 6) + 
    facet_grid(cluster ~ geneset) + 
    scale_color_viridis_c() + 
    theme_bw(24) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
    ggtitle(paste("Volcano plots by gene sets:", jmark))
})
print(m.vol.gsets)

gsets.filt <- c("Bcell", "Erythroblast", "HSCs", "Neutrophil")
m.vol.gsets.filt <- lapply(jmarks, function(jmark){
  m <- ggplot(fits.bygenesets.long %>% filter(abs(logLambda) < 10 & mark == jmark & geneset %in% gsets.filt), aes(x = logLambda, y = -log10(pval), color = logintercept)) +
    geom_point_rast(alpha = 0.5, size = 6) + 
    facet_grid(cluster ~ geneset) + 
    scale_color_viridis_c() + 
    theme_bw(24) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + 
    ggtitle(paste("Volcano plots by gene sets:", jmark))
})
print(m.vol.gsets.filt)

jgenes <- c("Hlf", "Hbb-bh2", "Pax5", "S100a7a")

for (jmark in jmarks){
  for (jgene in jgenes){
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

gdat.all <- lapply(jmarks, function(jmark){
  gene.mean <- apply(tss.mats.filt.fromref.cellfilt[[jmark]], 1, mean)
  gene.sd <- apply(tss.mats.filt.fromref.cellfilt[[jmark]], 1, sd)
  gene.cv <- gene.sd / gene.mean
  gdat <- data.frame(gene = rownames(tss.mats.filt.fromref.cellfilt[[jmark]]), 
                     gene.mean = gene.mean, 
                     gene.sd = gene.sd, 
                     gene.cv = gene.cv, 
                     mark = jmark, 
                     stringsAsFactors = FALSE)
}) %>%
  bind_rows()

head(gdat.all)

m1 <- ggplot(gdat.all, aes(x = gene.mean, y = gene.sd ^ 2)) + 
  geom_point_rast(alpha = 0.2) + 
  geom_density_2d() + 
  facet_wrap(~mark) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

m2 <- ggplot(gdat.all, aes(y = gene.cv ^ 2, x = gene.mean)) + 
  geom_abline(slope = -1, intercept = 0, linetype = "dotted", alpha = 0.3) + 
  scale_x_log10() + scale_y_log10() + 
  geom_point_rast(alpha = 0.2) + 
  facet_wrap(~mark) + 
  xlab("Mean Across Cells") + ylab("CV^2") + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

print(m1)
print(m2)


dev.off()

# if (!file.exists(outfits.wrangled)){
  save(jfits.long, fits.bygenesets.long, file = outfits.wrangled)
# }


# Show overdispersion ---------------------------------------------------------------

# jmark <- "H3K4me3"


# Show fits of specific genes ---------------------------------------------
# show fit of neutro-specific gene


# RefitPoissonForPlot <- function(jrow, cnames, datannots.filt.mark, ncuts.cell.mark){
#   jfit <- FitGlmRowClusters(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, returnobj = TRUE)
#   ci <- confint(jfit$fit.full)
#   
#   params.fc.dat <- data.frame(param = rownames(ci), logLambdaLower = ci[, 1], logLambdaUpper = ci[, 2]) %>%
#     filter(param != "(Intercept)") %>%
#     rowwise() %>%
#     mutate(Cluster = ifelse(param == "(Intercept)", "aHSPCs", gsub("Cluster", "", param)),
#            logLambda = mean((logLambdaLower + logLambdaUpper) / 2))
#   
#   params.int.dat <- data.frame(param = rownames(ci)[[1]], logLambdaLower = ci[1, 1], logLambdaUpper = ci[1, 2]) %>%
#     rowwise() %>%
#     mutate(Cluster = ifelse(param == "(Intercept)", "aHSPCs", gsub("Cluster", "", param)),
#            logLambda = mean((logLambdaLower + logLambdaUpper) / 2))
#   
#   params.mean.dat <- params.fc.dat %>%
#     mutate(logLambda = params.int.dat$logLambda + logLambda,
#            logLambdaLower = params.int.dat$logLambdaLower + logLambdaLower,
#            logLambdaUpper = params.int.dat$logLambdaUpper + logLambdaUpper) %>%
#     bind_rows(., params.int.dat)
#   
#   # plot raw 
#   input.dat <- jfit$dat.input %>% rowwise() %>% mutate(logLambda = log(ncuts) - log(ncuts.total))
#   return(list(input.dat = input.dat, params.mean.dat = params.mean.dat, params.int.dat = params.int.dat))
# }

# 
# subset(jfits.dat.lst$H3K4me3, grepl("Hbb", gene))
# 
# head(subset(fits.bygenesets.long, geneset == "Erythroblast" & cluster == "ClusterErythroblasts") %>% arrange(desc(logLambda)))
# head(subset(fits.bygenesets.long, geneset == "Bcell" & cluster == "ClusterBcells") %>% arrange(desc(logLambda)), 100)
# 
# # jgene <- "S100a7a"
# jgene <- "Hlf"
# # jgene <- "Sox6"
# # jgene <- "Hbb-bh2"
# # jgene <- "Add2"
# 
# jmark <- "H3K4me3"
# jgene <- "Hbb-bh2"
# 
# 

# Scrap -------------------------------------------------------------------

# 
# jfit <- FitGlmRowClusters(jrow, cnames, dat.annots.filt.mark, ncuts.cells.mark, returnobj = TRUE)
# 
# ci <- confint(jfit$fit.full)
# 
# # get fits
# 
# params.fc.dat <- data.frame(param = rownames(ci), logLambdaLower = ci[, 1], logLambdaUpper = ci[, 2]) %>%
#   filter(param != "(Intercept)") %>%
#   rowwise() %>%
#   mutate(Cluster = ifelse(param == "(Intercept)", "aHSPCs", gsub("Cluster", "", param)),
#          logLambda = mean((logLambdaLower + logLambdaUpper) / 2))
# 
# params.int.dat <- data.frame(param = rownames(ci)[[1]], logLambdaLower = ci[1, 1], logLambdaUpper = ci[1, 2]) %>%
#   rowwise() %>%
#   mutate(Cluster = ifelse(param == "(Intercept)", "aHSPCs", gsub("Cluster", "", param)),
#          logLambda = mean((logLambdaLower + logLambdaUpper) / 2))
# 
# params.mean.dat <- params.fc.dat %>%
#   mutate(logLambda = params.int.dat$logLambda + logLambda,
#          logLambdaLower = params.int.dat$logLambdaLower + logLambdaLower,
#          logLambdaUpper = params.int.dat$logLambdaUpper + logLambdaUpper) %>%
#   bind_rows(., params.int.dat)
# 
# # plot raw 
# input.dat <- jfit$dat.input %>% rowwise() %>% mutate(logLambda = log(ncuts) - log(ncuts.total))
# ggplot(input.dat, aes(x = Cluster, y = logLambda)) + 
#   geom_jitter(width = 0.1, height = 0, alpha = 0.1) + 
#   geom_errorbar(mapping = aes(ymin = logLambdaLower, ymax = logLambdaUpper), data = params.mean.dat, width = 0.1) + 
#   theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   xlab("" ) + ylab("log(ncuts) - log(total cuts)") + 
#   ggtitle(paste(jgene, jbin, jmark)) + 
#   geom_hline(yintercept = params.int.dat$logLambda, linetype = "dotted")
# 
# 
# 
# 

