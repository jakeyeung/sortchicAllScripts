# Jake Yeung
# Date of Creation: 2020-11-08
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/primetime_plots/2-ncells_per_well_K562_plots.R
# 



rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(ggrastr)


# Load genomeiwde summaries -----------------------------------------------

# inf.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins/GenomeWideFits.2020-07-22.rds"

jdate <- "2020-07-24"
infrds <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins/GenomeWideFits.", jdate, ".logncells.rds")
jsub.sum <- readRDS(infrds)

# jsub.sum <- readRDS(inf.meta)
dat.meta <- jsub.sum %>%
  rowwise() %>%
  mutate(ncells.lin = ncells,
         ncells = log(ncells))

outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/primetime2/K562_ncells_analysis.", Sys.Date(), ".pdf")

pdf(file = outpdf, useDingbats = FALSE)


# ggplot(dat.meta %>% filter(ncells %in% c(1,2)), aes(x = spikeinconcFactor, y = log2(chromocounts / spikeincounts), fill = as.character(ncells))) +
#   geom_boxplot() +
#   facet_grid(mark ~ conc) +
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# ggplot(dat.meta %>% filter(spikeinconc == 50000 & conc == "37U"), aes(x = log(ncells.lin), y = log2(chromocounts / spikeincounts))) +
#   geom_point() + 
#   facet_grid(mark ~ conc) +
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Get g ood cells ---------------------------------------------------------


mincounts <- 1000



# Redo fits ---------------------------------------------------------------

jmarks <- c("H3K4me3", "H3K27me3"); names(jmarks) <- jmarks

# good.concs <- c(350, 1500, 3000, 12000, 25000); names(good.concs) <- good.concs
good.concs <- sort(as.integer(unique(as.character(dat.meta$spikeinconcFactor)))); names(good.concs) <- good.concs
mnaseus <- c("37U", "75U"); names(mnaseus) <- mnaseus



jsub <- subset(dat.meta %>% filter(spikeinconcFactor %in% good.concs & mark %in% jmarks & ncells.lin > 0 & chromocounts > mincounts)) %>%
  rowwise() %>%
  mutate(ncells = as.integer(ncells.lin),
         logncells = log(ncells),
         logratio = log(chromocounts / spikeincounts),
         log2ncells = log2(ncells),
         log2ratio = log2(chromocounts / spikeincounts), 
         log2chromo = log2(chromocounts), 
         logchromo = log(chromocounts))

# remove an outlier at conc 6000, 37U H3K4me3
print(dim(jsub))
bad.cell <- subset(jsub, mark == "H3K4me3" & spikeinconc == 6000 & conc == "37U" & log2ratio > 7)$samp

jsub <- subset(jsub, samp != bad.cell)
# jsub <- subset(jsub, mark != "H3K4me3" & spikeinconc != 6000 & log2ratio < 5)
print(dim(jsub))


jfits.byconc <- lapply(jmarks, function(jmark){
  lapply(good.concs, function(jconc){
    lapply(mnaseus, function(mnaseu){
      print(paste(jmark, jconc, mnaseu))
      jsubsub <- jsub %>% filter(mark == jmark & spikeinconcFactor == jconc & conc == mnaseu)
      print(nrow(jsubsub))
      # m <- ggplot(jsubsub, aes(x = logncells, y = logratio)) + 
      #   geom_point() + ggtitle(paste(jmark, jconc, mnaseu))
      #   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      # print(m)
      jfit <- lm(formula = logratio ~ 1 + logncells, data = jsubsub) 
      jdat <- data.frame(param = rownames(summary(jfit)$coefficients), summary(jfit)$coefficients, stringsAsFactors = FALSE)
      jdat$mark <- jmark
      jdat$conc <- jconc
      jdat$mnaseu <- mnaseu
      jdat$N <- nrow(jsubsub)
      jdat$Response <- "logratio"
      return(jdat)
    }) %>%
      bind_rows()
  }) %>%
    bind_rows()
}) %>%
  bind_rows()

jfits.byconc.nospikeins <- lapply(jmarks, function(jmark){
  lapply(good.concs, function(jconc){
    lapply(mnaseus, function(mnaseu){
      print(paste(jmark, jconc, mnaseu))
      jsubsub <- jsub %>% filter(mark == jmark & spikeinconcFactor == jconc & conc == mnaseu)
      print(nrow(jsubsub))
      # m <- ggplot(jsubsub, aes(x = logncells, y = logratio)) + 
      #   geom_point() + ggtitle(paste(jmark, jconc, mnaseu))
      #   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      # print(m)
      jfit <- lm(formula = logchromo ~ 1 + logncells, data = jsubsub) 
      jdat <- data.frame(param = rownames(summary(jfit)$coefficients), summary(jfit)$coefficients, stringsAsFactors = FALSE)
      jdat$mark <- jmark
      jdat$conc <- jconc
      jdat$mnaseu <- mnaseu
      jdat$N <- nrow(jsubsub)
      jdat$Response <- "logchromo"
      return(jdat)
    }) %>%
      bind_rows()
  }) %>%
    bind_rows()
}) %>%
  bind_rows()

ggplot(subset(jfits.byconc, param == "logncells"), aes(x = as.character(conc), y = Estimate)) + 
  geom_boxplot() + 
  facet_grid(mark ~ mnaseu)

ggplot(subset(jfits.byconc, param == "logncells"), aes(x = Estimate)) + 
  geom_histogram() +
  facet_wrap(~mark) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(subset(jfits.byconc.nospikeins, param == "logncells"), aes(x = Estimate)) + 
  geom_histogram() +
  facet_wrap(~mark) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# plot(density(subset(jfits.byconc, param == "logncells")$Estimate))

# plot some bad fits
outliers <- subset(jfits.byconc, param == "logncells" & Estimate > 1.2) %>%
  arrange(desc(abs(Estimate)))

outliers <- subset(jfits.byconc.nospikeins, param == "logncells" & Estimate < 0.5) %>%
  arrange(desc(abs(Estimate)))


# Plot outliers -----------------------------------------------------------

i <- 1
jmark <- outliers$mark[[i]]
jconc <- outliers$conc[[i]]
jmnaseu <- outliers$mnaseu[[i]]
jslope <- signif(outliers$Estimate[[i]], digits = 3)

jsubsub.check <- jsub %>% filter(mark == jmark & spikeinconcFactor == jconc & conc == jmnaseu)

refit <- lm(formula = logratio ~ 1 + logncells, data = jsubsub.check)
jslope.new <- signif(refit$coefficients[[2]], digits = 3)

m <- ggplot(jsubsub.check, aes(x = logncells, y = logratio)) + geom_point() + ggtitle(paste(jmark, jconc, jmnaseu, jslope, jslope.new)) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_smooth(se = FALSE, method = "lm")
print(m)

m <- ggplot(jsubsub.check, aes(x = logncells, y = logchromo)) + geom_point() + ggtitle(paste(jmark, jconc, jmnaseu, jslope, jslope.new)) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_smooth(se = FALSE, method = "lm")
print(m)



# Fit across all concentrations and mnaseUs -------------------------------


fits.all.spikeins <- lapply(jmarks, function(jmark){
  jsuball <- subset(jsub, mark == jmark)
  fitall <- lm(formula = logratio ~ 1 + spikeinconcFactor + conc + logncells, data = jsuball)
  datall <- data.frame(param = rownames(summary(fitall)$coefficients), summary(fitall)$coefficients, stringsAsFactors = FALSE) 
  datall$mark <- jmark
  return(datall)
}) %>%
  bind_rows()


fits.all.nospikeins <- lapply(jmarks, function(jmark){
  jsuball <- subset(jsub, mark == jmark)
  fitall <- lm(formula = logchromo ~ 1 + spikeinconcFactor + conc + logncells, data = jsuball)
  datall <- data.frame(param = rownames(summary(fitall)$coefficients), summary(fitall)$coefficients, stringsAsFactors = FALSE) 
  datall$mark <- jmark
  return(datall)
}) %>%
  bind_rows()


# Show all slopes ---------------------------------------------------------

m.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(subset(jsub, mark == jmark), aes(x = log2ncells, y = log2ratio)) + 
    geom_point() + 
    theme_bw() + 
    facet_grid(conc ~ spikeinconcFactor) + 
    ggtitle(jmark) + 
    geom_smooth(method = "lm", se = TRUE) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(breaks = log2(seq(5)), labels = seq(5)) + 
    xlab("Number of Cells") + ylab("log2(cuts/spikeins)")
  return(m)
})
print(m.lst)

JFuncs::multiplot(m.lst$H3K4me3, m.lst$H3K27me3, cols = 1)


m.nospikeins.lst <- lapply(jmarks, function(jmark){
  m <- ggplot(subset(jsub, mark == jmark), aes(x = log2ncells, y = log2chromo)) + 
    geom_point() + 
    theme_bw() + 
    facet_grid(conc ~ spikeinconcFactor) + 
    ggtitle(jmark) + 
    geom_smooth(method = "lm", se = TRUE) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(breaks = log2(seq(5)), labels = seq(5)) + 
    xlab("Number of Cells") + ylab("log2(cuts)")
  return(m)
})
print(m.nospikeins.lst)

JFuncs::multiplot(m.nospikeins.lst$H3K4me3, m.nospikeins.lst$H3K27me3, cols = 1)




# Summarize the fits ------------------------------------------------------

ggplot(subset(jfits.byconc, param == "logncells"), aes(x = Estimate, fill = mark)) + 
  geom_histogram(alpha = 0.25) + 
  coord_cartesian(xlim = c(0, 2)) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(subset(jfits.byconc.nospikeins, param == "logncells"), aes(x = Estimate, fill = mark)) + 
  geom_histogram(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

jfits.byconc.merge <- rbind(jfits.byconc, jfits.byconc.nospikeins)

jtmp <- subset(jfits.byconc.merge, param == "logncells") %>% mutate(concmnase = interaction(conc, mnaseu))

jfits.byconc.merge.wide <- reshape2::dcast(jtmp, concmnase ~ mark + Response, value.var = "Estimate")

ggplot(jfits.byconc.merge.wide, aes(x = H3K27me3_logratio, y = H3K27me3_logchromo)) + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jfits.byconc.merge.wide, aes(x = H3K4me3_logratio, y = H3K4me3_logchromo)) + 
  geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

library(ggrepel)
ggplot(jfits.byconc.merge.wide, aes(x = sqrt(H3K4me3_logratio * H3K27me3_logratio), y = sqrt(H3K4me3_logchromo + H3K27me3_logratio), label = concmnase)) + 
  geom_point()  + 
  geom_text_repel() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(subset(jfits.byconc.merge, param == "logncells"), aes(x = Estimate, fill = Response)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Slope")

for (jres in unique(jfits.byconc.merge$Response)){
  m <- ggplot(subset(jfits.byconc.merge, param == "logncells" & Response == jres), aes(x = Estimate, fill = Response)) + 
    geom_density(alpha = 0.25) + 
    ggtitle(jres) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    xlab("Slope") + 
    geom_vline(xintercept = 1, linetype = "dotted")
    xlim(c(0, 2))
  print(m)
}


# Do a fit for K4me3 ------------------------------------------------------

jconc <- 12000
jmnaseu <- "75U"

m.lst <- lapply(jmarks, function(jmark){
  jsubsub.check <- jsub %>% filter(mark == jmark & spikeinconcFactor == jconc & conc == jmnaseu)
  jfits.check <- subset(jfits.byconc.merge, param == "logncells" & conc == jconc & mnaseu == jmnaseu & mark == jmark)
  jslope.ratio <- signif(subset(jfits.check, Response == "logratio")$Estimate[[1]], digits = 3)
  jslope.chromo <- signif(subset(jfits.check, Response == "logchromo")$Estimate[[1]], digits = 3)
  m1 <- ggplot(jsubsub.check, aes(x = log2ncells, y = log2ratio)) + geom_point() + ggtitle(paste(jmark, jconc, jmnaseu, jslope.ratio)) + 
    theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_smooth(se = TRUE, method = "lm") + 
    xlab("Number of Cells") + ylab("log2(cuts/spikeins)") + 
    scale_x_continuous(breaks = log2(seq(5)), labels = seq(5)) 

  m2 <- ggplot(jsubsub.check, aes(x = log2ncells, y = log2chromo)) + geom_point() + ggtitle(paste(jmark, jconc, jmnaseu, jslope.chromo)) + 
    theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_smooth(se = TRUE, method = "lm") + 
    xlab("Number of Cells") + ylab("log2(cuts)") + 
    scale_x_continuous(breaks = log2(seq(5)), labels = seq(5)) 
  return(list(m1, m2))
})

print(m.lst)


# Do a fit for K27me3 ------------------------------------------------------



# Calculate variance ------------------------------------------------------

jsub.sd <- jsub %>%
  group_by(mark, conc, spikeinconcFactor) %>%
  summarise(log2sdchromo = sd(log2chromo),
            log2sdratio = sd(log2ratio),
            cvchromo = sd(chromocounts) / mean(chromocounts),
            cvspikein = sd(chromocounts / spikeincounts) / mean(chromocounts / spikeincounts))

jsub.sd.long <- jsub.sd %>%
  melt()

colnames(jsub.sd.long) <- c("mark", "conc", "spikeinconcFactor", "Response", "StdDev")

ggplot(jsub.sd.long %>% filter(!grepl("^cv", Response)), aes(x = StdDev, fill = Response)) + 
  geom_density(alpha = 0.25) +
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  xlab("StdDev")

ggplot(jsub.sd.long %>% filter(!grepl("^cv", Response)), aes(x = Response, y = StdDev)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.1) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ylab("StdDev") + xlab("")

ggplot(jsub.sd.long %>% filter(!grepl("^cv", Response)), aes(x = forcats::fct_reorder(.f = Response, .x = StdDev, .fun = mean, .desc = FALSE), y = StdDev)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.1) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ylab("StdDev") + xlab("")

ggplot(jsub.sd.long %>% filter(!grepl("^cv", Response)), aes(x = forcats::fct_reorder(.f = Response, .x = StdDev, .fun = mean, .desc = FALSE), y = StdDev)) + 
  geom_boxplot() + 
  facet_wrap(~mark) + 
  geom_jitter(width = 0.1) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ylab("StdDev") + xlab("")


ggplot(jsub.sd.long %>% filter(grepl("^cv", Response)), aes(x = Response, y = StdDev)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.1) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ylab("CV") + xlab("")

ggplot(jsub.sd.long %>% filter(grepl("^cv", Response)), aes(x = forcats::fct_reorder(.f = Response, .x = StdDev, .fun = mean, .desc = FALSE), y = StdDev)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.1) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ylab("CV") + xlab("")

ggplot(jsub.sd.long %>% filter(grepl("^cv", Response)), aes(x = forcats::fct_reorder(.f = Response, .x = StdDev, .fun = mean, .desc = FALSE), y = StdDev)) + 
  geom_boxplot() + 
  facet_wrap(~mark) + 
  geom_jitter(width = 0.1) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ylab("CV") + xlab("")


print("Done")

# do statistical test 
# p-value < 
jrow <- jsub.sd[1, ]
jcheck <- subset(jsub, mark == "H3K4me3" & conc == "37U" & 350)
jout <- pf( var(jsub$log2chromo) / var(jsub$log2ratio), df1 = nrow(jsub) - 1, df2 = nrow(jsub) - 1, lower.tail = TRUE)
var.test(jsub$log2chromo, jsub$log2ratio, ratio = 1, alternative = "two.sided")

dev.off()


# Check outlier -----------------------------------------------------------

# H3K4me3 has high CV estimates for 2 fits, why? remove them if needed 


x <- jsub.sd.long %>% filter(grepl("^cv", Response))

ggplot(jsub.sd.long %>% filter(grepl("^cv", Response)), aes(x = forcats::fct_reorder(.f = Response, .x = StdDev, .fun = mean, .desc = FALSE), y = StdDev)) + 
  geom_boxplot(outlier.shape = NA) + 
  facet_wrap(~mark) + 
  geom_jitter(width = 0.1) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ylab("CV") + xlab("")



# 
# 
# # Fit full model  ---------------------------------------------------------
# 
# fit.spikeins.full <- lm(formula = logratio ~ 1 + spikeinconcFactor + mark + logncells:mark, data = jsub)
# fit.nospikeins.full <- lm(formula = logchromo ~ 1 + spikeinconcFactor + mark + logncells:mark, data = jsub)
# 
# # calculate variance
# 
# # do test 
# jtest.var <- var.test(fit.spikeins.full, fit.nospikeins.full, ratio = 1)
# 
# # confirm manually 
# var1 <- sum(fit.spikeins.full$residuals ^ 2) / nrow(jsub)
# var2 <- sum(fit.nospikeins.full$residuals ^ 2) / nrow(jsub)
# jout <- 1 - pf( var1 / var2, df1 = fit.spikeins.full$df.residual, df2 = fit.nospikeins.full$df.residual, lower.tail = TRUE)
# 
# 
# fits.all.spikeins <- lapply(jmarks, function(jmark){
#   jsuball <- subset(jsub, mark == jmark)
#   fitall <- lm(formula = logratio ~ 1 + spikeinconcFactor + conc + logncells, data = jsuball)
#   datall <- data.frame(param = rownames(summary(fitall)$coefficients), summary(fitall)$coefficients, stringsAsFactors = FALSE) 
#   datall$mark <- jmark
#   return(datall)
# }) %>%
# bind_rows()
# 
# 
# fits.all.nospikeins <- lapply(jmarks, function(jmark){
#   jsuball <- subset(jsub, mark == jmark)
#   fitall <- lm(formula = logchromo ~ 1 + spikeinconcFactor + conc + logncells, data = jsuball)
#   datall <- data.frame(param = rownames(summary(fitall)$coefficients), summary(fitall)$coefficients, stringsAsFactors = FALSE) 
#   datall$mark <- jmark
#   return(datall)
# }) %>%
#   bind_rows()
