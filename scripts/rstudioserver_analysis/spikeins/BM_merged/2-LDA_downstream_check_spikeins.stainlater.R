# Jake Yeung
# Date of Creation: 2020-10-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged/2-LDA_downstream_check_spikeins.R
# 



rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)
library(scchicFuncs)

library(topicmodels)

library(ggrepel)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)


jmarks <- c("H3K27me3"); names(jmarks) <- jmarks


# Plot chromo spikeins across plates k27me3 -------------------------------

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_mouse.BMround2all.stainlater.varfilt"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_mouse.BMround2all.stainlater.varfilt"
inrdata <- file.path(indir, "LDA_downstream_objects.RData")

load(inrdata, v=T)

dat.merge$stype <- factor(as.character(dat.merge$stype), levels = c("LSK", "LinNeg", "Unenriched"))
dat.merge$mark <- jmarks[[1]]

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.merge, aes(x = log2(chromocounts / spikeincounts), fill = stype)) + 
  geom_density(alpha = 0.25)  + 
  theme_bw() + 
  scale_fill_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark)

ggplot(dat.merge %>% filter(mark == "H3K27me3"), aes(x = log2(chromocounts / spikeincounts), fill = stype)) + 
  geom_density(alpha = 0.25)  + 
  theme_bw() + 
  scale_fill_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~plate)


outf <- file.path(outdir, "spikein_stats_across_stype.pdf")
# assertthat::assert_that(!file.exists(outf))

pdf(file = outf, useDingbats = FALSE)
for (jmark in jmarks){
  m <- ggplot(dat.merge %>% filter(mark == jmark), aes(y = log2(chromocounts / spikeincounts), x = stype)) + 
    geom_boxplot() + 
    theme_bw() + 
    scale_fill_manual(values = cbPalette) + 
    facet_wrap(~experi) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), strip.text.x = element_text(size = 5)) + 
    ggtitle(jmark)
  print(m)
}


for (jmark in jmarks){
  m <- ggplot(dat.merge %>% filter(mark == jmark), aes(y = log2(chromocounts), x = stype)) + 
    geom_boxplot() + 
    theme_bw() + 
    scale_fill_manual(values = cbPalette) + 
    facet_wrap(~experi) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), strip.text.x = element_text(size = 5)) + 
    ggtitle(jmark)
  print(m)
}


for (jmark in jmarks){
  m <- ggplot(dat.merge %>% filter(mark == jmark), aes(y = cell.var.within.sum.norm, x = stype)) + 
    geom_boxplot() + 
    theme_bw() + 
    scale_fill_manual(values = cbPalette) + 
    facet_wrap(~experi) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), strip.text.x = element_text(size = 5)) + 
    ggtitle(jmark)
  print(m)
}

dev.off()

# plot UMAPs by plate -----------------------------------------------------

outf2 <- file.path(outdir, "spikein_stats_across_stype_umaps.pdf")
pdf(file = outf2, useDingbats = FALSE)
for (jmark in jmarks){
  jexperis <- unique(subset(dat.merge, mark == jmark)$experi)
  for (jexperi in jexperis){
    m <- ggplot(dat.merge %>% filter(mark == jmark & experi == jexperi), aes(x = umap1, y = umap2, color = log2(chromocounts / spikeincounts))) + 
      geom_point() + 
      theme_bw() + 
      scale_color_viridis_c() + 
      facet_wrap(~experi) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), strip.text.x = element_text(size = 5)) + 
      ggtitle(jmark, jexperi)
    print(m)
  }
}
dev.off()



# Fit models --------------------------------------------------------------

DoFits <- function(jsub){
  jsub$stype <- factor(as.character(jsub$stype), levels = c("Unenriched", "LinNeg", "LSK"))
  
  jfit <- lm(formula = l2r ~ 1 + experi + stype, data = jsub)
  jfit.null <- lm(formula = l2r ~ 1 + experi, data = jsub)
  jcompare <- anova(jfit.null, jfit)
  jfit.effects <- data.frame(param = rownames(summary(jfit)$coefficients), summary(jfit)$coefficients, stringsAsFactors = FALSE) %>%
    filter(!grepl("^experi|Intercept", param))
  jfit.int <- data.frame(param = rownames(summary(jfit)$coefficients), summary(jfit)$coefficients, stringsAsFactors = FALSE) %>%
    filter(param == "(Intercept)") %>%
    rowwise() 
  jint <- jfit.int$Estimate
  jint.se <- jfit.int$Std..Error
  
  jfit.int <- data.frame(param = rownames(summary(jfit)$coefficients), summary(jfit)$coefficients, stringsAsFactors = FALSE) %>%
    filter(param == "(Intercept)") %>%
    rowwise() %>%
    mutate(est = jint,
           est.se = jint.se)
  
  jfit.dat <- jfit.effects %>%
    rowwise() %>%
    mutate(est = Estimate + jint,
           est.se = sqrt(Std..Error ^ 2 + jint.se ^ 2))
  
  jfit.merge <- rbind(jfit.dat, jfit.int)
  jfit.merge$param <- gsub("stype", "", jfit.merge$param)
  jfit.merge$param <- gsub("(Intercept)", "Unenriched", jfit.merge$param)
  return(list(jfit = jfit, jfit.null = jfit.null, jcompare = jcompare, jfit.merge = jfit.merge))
  
}


jouts.l2r <- lapply(jmarks, function(jmark){
  print(jmark)
  jsub <- subset(dat.merge, mark == jmark) %>%
    rowwise() %>%
    mutate(l2r = log2(chromocounts / spikeincounts))
  jout <- DoFits(jsub)
  jout$jfit.merge$mark <- jmark
  return(jout)
}) 

jouts.cuts <- lapply(jmarks, function(jmark){
  jsub <- subset(dat.merge, mark == jmark) %>%
    rowwise() %>%
    mutate(l2r = log2(chromocounts))
  jout <- DoFits(jsub)
  jout$jfit.merge$mark <- jmark
  return(jout)
}) 

jouts.var <- lapply(jmarks, function(jmark){
  jsub <- subset(dat.merge, mark == jmark) %>%
    rowwise() %>%
    mutate(l2r = cell.var.within.sum.norm)
  jout <- DoFits(jsub)
  jout$jfit.merge$mark <- jmark
  return(jout)
}) 


jfits.merge.l2r <- lapply(jouts.l2r, function(x) x$jfit.merge) %>%
  bind_rows()

jfits.merge.cuts <- lapply(jouts.cuts, function(x) x$jfit.merge) %>%
  bind_rows()

jfits.merge.var <- lapply(jouts.var, function(x) x$jfit.merge) %>%
  bind_rows()

outf3 <- file.path(outdir, "spikein_stats_across_stype_estimates.pdf")
pdf(file = outf3, useDingbats = FALSE)
  m <- ggplot(jfits.merge.l2r, aes(x = param, y = est, ymin = est - 1.96 * est.se, ymax = est + 1.96 * est.se)) + 
    geom_point() + 
    geom_errorbar() + 
    facet_wrap(~mark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    xlab("") + ylab("Estiamtes for log2(cuts / spikeins)")
  print(m)
  m <- ggplot(jfits.merge.cuts, aes(x = param, y = est, ymin = est - 1.96 * est.se, ymax = est + 1.96 * est.se)) + 
    geom_point() + 
    geom_errorbar() + 
    facet_wrap(~mark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    xlab("") + ylab("Estimates for log2(cuts)")
  print(m)
  m <- ggplot(jfits.merge.var, aes(x = param, y = est, ymin = est - 1.96 * est.se, ymax = est + 1.96 * est.se)) + 
    geom_point() + 
    geom_errorbar() + 
    facet_wrap(~mark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    xlab("") + ylab("Estimates for var")
  print(m)
dev.off()

# Plot fits  --------------------------------------------------------------



# 
# for (jmark in jmarks){
#   m <- ggplot(dat.merge %>% filter(mark == jmark), aes(y = log2(chromocounts), x = stype)) + 
#     geom_boxplot() + 
#     theme_bw() + 
#     scale_fill_manual(values = cbPalette) + 
#     facet_wrap(~experi) + 
#     theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), strip.text.x = element_text(size = 5)) + 
#     ggtitle(jmark)
#   print(m)
# }
# 
# 
# for (jmark in jmarks){
#   m <- ggplot(dat.merge %>% filter(mark == jmark), aes(y = cell.var.within.sum.norm, x = stype)) + 
#     geom_boxplot() + 
#     theme_bw() + 
#     scale_fill_manual(values = cbPalette) + 
#     facet_wrap(~experi) + 
#     theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), strip.text.x = element_text(size = 5)) + 
#     ggtitle(jmark)
#   print(m)
# }




# Fit model  --------------------------------------------------------------


