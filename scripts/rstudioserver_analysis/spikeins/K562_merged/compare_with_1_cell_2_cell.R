# Jake Yeung
# Date of Creation: 2020-09-17
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_merged/compare_with_1_cell_2_cell.R
# Calculate effect sizes of 1 and 2 cell experiment 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_K562.round2.AllMerged.WithDoubleCells/compare_merged_with_1_2_cells.pdf")

pdf(outpdf, useDingbats = FALSE)

# Load different spikeins measure counts in each chromosome ------------------------------------------------

spikeinchromo <- "J02459.1"
jchromos <- paste("", seq(19), sep = "")

hubpath <- "/home/jyeung/hub_oudenaarden"

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/spikein/fastqs/tagged_bams.scmo3.contigfixed/countTablesAndRZr1only_ByChromo"


infs.lst <- list.files(indir, pattern = "*.csv", full.names = TRUE)
names(infs.lst) <- sapply(infs.lst, basename)


jmark <- "H3K4me3"
jconc <- "37U"


jmarks.lst <- c("H3K4me3", "H3K27me3"); names(jmarks.lst) <- jmarks.lst
jconcs.lst <- c("37U", "75U"); names(jconcs.lst) <- jconcs.lst

inf <- file.path(hubpath, paste0("jyeung/data/scChiC/spikein/fastqs/tagged_bams.scmo3.contigfixed/countTablesAndRZr1only_ByChromo/PZ-K562-", jmark, "-spikein-", jconc, ".scmo3.again_contigfixed.tagged.countTable.ByChromo.csv"))

assertthat::assert_that(file.exists(inf))


# Load object --------------------------------------------------------------

jdate <- "2020-07-24"
infrds <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins/GenomeWideFits.", jdate, ".logncells.rds")

dat.sum <- readRDS(infrds)

ggplot(dat.sum %>% filter(ncells %in% c(1,2)), aes(x = spikeinconcFactor, y = log2(chromocounts / spikeincounts), fill = as.character(ncells))) + 
  geom_boxplot() + 
  # geom_point()  + 
  facet_grid(mark ~ conc) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# pick some good concentrations -------------------------------------------

jmarks <- c("H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

good.concs <- c(350, 1500, 3000, 12000, 25000); names(good.concs)

jsub <- subset(dat.sum %>% filter(ncells %in% c(1, 2) & spikeinconcFactor %in% good.concs & mark %in% jmarks)) %>%
  ungroup() %>%
  mutate(ncells = as.character(ncells))

ggplot(jsub, aes(x = spikeinconcFactor, y = log2(chromocounts / spikeincounts), fill = ncells)) + 
  geom_boxplot() + 
  # geom_point()  + 
  facet_grid(mark ~ conc) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# estimate effect size 

jfits <- lapply(jmarks, function(jmark){
  jfit <- lm(formula = log2(chromocounts / spikeincounts) ~ spikeinconcFactor + ncells, data = jsub %>% filter(mark == jmark)) 
  jdat <- data.frame(param = rownames(summary(jfit)$coefficients), summary(jfit)$coefficients, stringsAsFactors = FALSE)
  jdat$mark <- jmark
  return(jdat)
}) %>%
  bind_rows()

ggplot(jfits %>% filter(param == "ncells2"), aes(x = mark, y = Estimate, ymin = Estimate - 1.96 * Std..Error, ymax = Estimate + 1.96 * Std..Error)) + 
  geom_point() + 
  geom_errorbar() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# fit factor b yfactor  

jfits.byconc <- lapply(jmarks, function(jmark){
  lapply(good.concs, function(jconc){
    jfit <- lm(formula = log2(chromocounts / spikeincounts) ~ 1 + ncells, data = jsub %>% filter(mark == jmark & spikeinconcFactor == jconc)) 
    jdat <- data.frame(param = rownames(summary(jfit)$coefficients), summary(jfit)$coefficients, stringsAsFactors = FALSE)
    jdat$mark <- jmark
    jdat$conc <- jconc
    return(jdat)
  }) %>%
    bind_rows()
}) %>%
  bind_rows()

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
ggplot(jfits.byconc %>% filter(param == "ncells2"), 
       aes(x = mark, y = Estimate, ymin = Estimate - 1.96 * Std..Error, ymax = Estimate + 1.96 * Std..Error)) + 
  geom_boxplot() + 
  geom_jitter(mapping = aes(color = as.factor(conc)), width = 0.125) + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


 
dev.off()


# Dat filt  ---------------------------------------------------------------
# 
# 
# dat.filt.long <- lapply(jmarks.lst, function(jmark){
#   print(jmark)
#   lapply(jconcs.lst, function(jconc){
#     print(jconc)
#     inf <- file.path(hubpath, paste0("jyeung/data/scChiC/spikein/fastqs/tagged_bams.scmo3.contigfixed/countTablesAndRZr1only_ByChromo/PZ-K562-", jmark, "-spikein-", jconc, ".scmo3.again_contigfixed.tagged.countTable.ByChromo.csv"))
#     print(inf)
#     dat.filt.long <- GetChromoCounts(inf)
#     dat.filt.long$mark <- jmark
#     dat.filt.long$conc <- jconc
#     return(dat.filt.long)
#   })  %>%
#     bind_rows()
# }) %>%
#   bind_rows()
# 
# dat.rz <- lapply(jmarks.lst, function(jmark){
#   print(jmark)
#   lapply(jconcs.lst, function(jconc){
#     print(jconc)
#     inf.rz <- file.path(hubpath, paste0("jyeung/data/scChiC/spikein/fastqs/tagged_bams.scmo3.contigfixed//RZcounts/PZ-K562-", jmark, "-spikein-", jconc, ".scmo3.again_contigfixed.tagged.LH_counts.demuxbugfixed_mergeplates.csv"))
#     print(inf.rz)
#     dat.filt.long <- ReadLH.SummarizeTA(inf.rz)
#     dat.filt.long$mark <- jmark
#     dat.filt.long$conc <- jconc
#     return(dat.filt.long)
#   })  %>%
#     bind_rows()
# }) %>%
#   bind_rows()
# 
# 
# 
# 
# # Get global changes ------------------------------------------------------
# 
# 
# 
# 
# 
# # Calculate fold changes global -------------------------------------------
# 
# # fit linear model
# 
# print(unique(jsub.sum$spikeinconc))
# 
# input.dat <- subset(jsub.sum, spikeinconc == 50000 & conc == "37U" & mark == "H3K4me3" & ncells > 0)
# 
# 
# 
# spikeinconc.vec <- sort(unique(jsub.sum$spikeinconc))
# conc.vec <- c("37U", "75U")
# jmarks.vec <- c("H3K4me3", "H3K27me3")
# 
# for (jspikeinconc in spikeinconc.vec){
#   for (jconc in conc.vec){
#     for (jmark in jmarks.vec){
#       input.dat <- subset(jsub.sum, spikeinconc == jspikeinconc & conc == jconc & mark == jmark & ncells > 0) %>%
#         rowwise() %>%
#         mutate(ncells.lin = ncells,
#                ncells = log(ncells))  # makes fits linear
#       # plot fits
#       f1 <- FitNormCountsToNcells.lm(input.dat, return.fit.obj = TRUE)
#       f2 <- FitNormCountsToNcells.glm(input.dat, return.fit.obj = TRUE)
#       # f2.ci <- confint.default(f2)
#       # jslope.glm.ln <- coefficients(f2)[["ncells"]]
#       jpred1 <- data.frame(ypred = predict(f1, input.dat, se.fit = FALSE), ncells = input.dat$ncells, chromocounts = input.dat$chromocounts, spikeinconc = input.dat$spikeinconc, spikeincounts = input.dat$spikeincounts, stringsAsFactors = FALSE)
#       jpred2 <- data.frame(ypred.unadj = predict(f2, input.dat, se.fit = FALSE), ncells = input.dat$ncells, chromocounts = input.dat$chromocounts, spikeinconc = input.dat$spikeinconc, spikeincounts = input.dat$spikeincounts, stringsAsFactors = FALSE) %>%
#         mutate(ypred = ypred.unadj - log(spikeincounts))
#       m.check1 <- ggplot(input.dat, 
#                          aes(x = ncells, y = log(chromocounts / spikeincounts))) + geom_point()  + 
#         theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#         geom_line(data = jpred1, mapping = aes(x = ncells, y = ypred)) + 
#         ggtitle(paste("Lm fit: SpikeInMole:", jspikeinconc, jconc, jmark)) + 
#         xlab("log(ncells)")
#       print(m.check1)
#       m.check2 <- ggplot(input.dat, 
#                          aes(x = ncells, y = log(chromocounts / spikeincounts))) + geom_point()  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#         geom_line(data = jpred2, mapping = aes(x = ncells, y = ypred)) + 
#         ggtitle(paste("GLM fit: SpikeInMole:", jspikeinconc, jconc, jmark)) + 
#         xlab("log(ncells)")
#       print(m.check2)
#       # plot linear
#       
#       # jfit.l1 <- lm(formula = chromocounts ~ ncells.lin, data = input.dat)
#       # jfit.l2 <- lm(formula = chromocounts / spikeincounts ~ ncells.lin, data = input.dat)
#       m.l1 <- ggplot(input.dat, aes(x = ncells.lin, y = chromocounts)) + 
#         geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#         ggtitle(paste("Linear fit on linear scale SpikeInMole:", jspikeinconc, jconc, jmark)) + 
#         geom_smooth(method = "lm", se = FALSE) 
#       m.l2 <- ggplot(input.dat, aes(x = ncells.lin, y = chromocounts / spikeincounts)) + 
#         geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#         ggtitle(paste("Linear fit on linear scale SpikeInMole:", jspikeinconc, jconc, jmark)) + 
#         geom_smooth(method = "lm", se = FALSE) 
#     }
#   }
# }
