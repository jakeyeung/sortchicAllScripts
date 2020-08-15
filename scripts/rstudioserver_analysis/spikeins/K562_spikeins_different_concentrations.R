# Jake Yeung
# Date of Creation: 2020-07-20
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_spikeins_different_concentrations.R
# 

rm(list=ls())

library(tidyverse)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)


# Functions ---------------------------------------------------------------

outrds <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins/GenomeWideFits.", Sys.Date(), ".logncells.rds")
outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins/GenomeWideFits.", Sys.Date(), ".logncells.pdf")


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

# Load data  --------------------------------------------------------------



dat.filt.long <- lapply(jmarks.lst, function(jmark){
  print(jmark)
  lapply(jconcs.lst, function(jconc){
    print(jconc)
    inf <- file.path(hubpath, paste0("jyeung/data/scChiC/spikein/fastqs/tagged_bams.scmo3.contigfixed/countTablesAndRZr1only_ByChromo/PZ-K562-", jmark, "-spikein-", jconc, ".scmo3.again_contigfixed.tagged.countTable.ByChromo.csv"))
    print(inf)
    dat.filt.long <- GetChromoCounts(inf)
    dat.filt.long$mark <- jmark
    dat.filt.long$conc <- jconc
    return(dat.filt.long)
  })  %>%
    bind_rows()
}) %>%
  bind_rows()

dat.rz <- lapply(jmarks.lst, function(jmark){
  print(jmark)
  lapply(jconcs.lst, function(jconc){
    print(jconc)
    inf.rz <- file.path(hubpath, paste0("jyeung/data/scChiC/spikein/fastqs/tagged_bams.scmo3.contigfixed//RZcounts/PZ-K562-", jmark, "-spikein-", jconc, ".scmo3.again_contigfixed.tagged.LH_counts.demuxbugfixed_mergeplates.csv"))
    print(inf.rz)
    dat.filt.long <- ReadLH.SummarizeTA(inf.rz)
    dat.filt.long$mark <- jmark
    dat.filt.long$conc <- jconc
    return(dat.filt.long)
  })  %>%
    bind_rows()
}) %>%
  bind_rows()


# pdfout <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins/initial_explore.", Sys.Date(), ".pdf")
# pdf(, useDingbats = FALSE)
pdf(outpdf, width = 1020/72, height = 815/72, useDingbats = FALSE)


# Remove low counts? ------------------------------------------------------

mincounts <- 1000
cellsfilt <- subset(dat.filt.long, mincounts >= mincounts)$samp
jtitle <- paste("Cells with Autosome Counts > 1000")

jsub <- dat.filt.long %>% filter(samp %in% cellsfilt)

# get ratio of chromocounts / spikeincounts


ggplot(jsub, aes(x = chromo, y = counts)) + geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1)  + 
  theme_bw() +  
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  scale_y_log10() + facet_wrap(~experi, nrow = 1) + 
  ggtitle(jtitle)


ggplot(jsub %>% filter(chromo %in% spikeinchromo), aes(x = chromo, y = counts, fill = experi)) + geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1)  + 
  theme_bw() +  
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  scale_y_log10() + 
  ggtitle(jtitle)


ggplot(jsub, aes(x = chromo, y = counts / totalcounts)) + geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  scale_y_log10() + 
  facet_grid(mark ~ conc) + 
  # facet_wrap(~experi, nrow = 1) + 
  xlab("") + 
  ylab("Counts / TotalCounts") + 
  ggtitle(jtitle)

ggplot(jsub, aes(x = chromo, y = counts / totalcounts, fill = experi)) + geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1)  + 
  theme_bw() + 
  theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom") + 
  scale_y_log10() + 
  xlab("") + 
  ylab("Counts / TotalCounts") + 
  ggtitle(jtitle)


ggplot(jsub, aes(x = chromo, y = counts / chromocounts)) + geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1)  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  scale_y_log10() + 
  facet_grid(mark ~ conc) + 
  # facet_wrap(~experi, nrow = 1) + 
  xlab("") + 
  ylab("Counts / AutosomeCounts") +
  ggtitle(jtitle) 


# Normalize by spikeins see if we get "double"  ---------------------------

ggplot(jsub, aes(x = chromo, y = log2(counts / spikeincounts), fill = experi)) + geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1)  + 
  theme_bw() + 
  theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom") + 
  xlab("") + 
  ylab("log2(Counts / SpikeIn)")

ggplot(jsub, aes(x = chromo, y = log2(counts / totalcounts), fill = experi)) + geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1)  + 
  theme_bw() + 
  theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom") + 
  xlab("") + 
  ylab("log2(Counts / TotalCounts)") + 
  ggtitle(jtitle)

ggplot(jsub, aes(x = chromo, y = log2(counts / chromocounts), fill = experi)) + geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.1)  + 
  theme_bw() + 
  theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom") + 
  xlab("") + 
  ylab("log2(Counts / AutosomeCounts)") + 
  ggtitle(jtitle)

# plot the density plots genomewide? 

ggplot(jsub %>% filter(chromo != spikeinchromo), aes(x = log2(counts / spikeincounts), fill = experi)) + 
  geom_density(alpha = 0.25)  + 
  theme_bw() + 
  theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom") + 
  xlab("") + 
  ylab("log2(Counts / SpikeIn)") + 
  facet_wrap(~chromo) + 
  ggtitle(jtitle)

ggplot(jsub %>% filter(chromo %in% jchromos), aes(x = log2(counts / chromocounts), fill = experi)) + 
  geom_density(alpha = 0.25)  + 
  theme_bw() + 
  theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), legend.position = "bottom") + 
  xlab("") + 
  ylab("log2(Counts / AutosomeCounts)") + 
  facet_wrap(~chromo) + 
  ggtitle(jtitle)


# Plot grid ---------------------------------------------------------------

dat.filt.long <- dat.filt.long %>%
  rowwise() %>%
  mutate(cellid = paste("cell", strsplit(samp, "_")[[1]][[2]], sep = ""),
         indx = strsplit(samp, "_")[[1]][[2]], 
         rowcoord = GetPlateCoord(cell = cellid, platecols = 24, is.zero.base = FALSE)[[1]],
         colcoord = GetPlateCoord(cell = cellid, platecols = 24, is.zero.base = FALSE)[[2]])

jexp <- "PZ-K562-H3K4me3-spikein-37U"

jsub.sum <- subset(dat.filt.long) %>%
  group_by(samp, experi, mark, conc, cellid, indx, rowcoord, colcoord) %>%
  summarise(totalcounts = unique(totalcounts),
            spikeincounts = unique(spikeincounts),
            chromocounts = unique(chromocounts))
  

jsub.exp <- subset(dat.filt.long, experi == jexp) %>%
  group_by(samp, experi, mark, conc, cellid, indx, rowcoord, colcoord) %>%
  summarise(totalcounts = unique(totalcounts),
         spikeincounts = unique(spikeincounts),
         chromocounts = unique(chromocounts))


# jsub.exp <- jsub.sum


# check wellname is correct
ggplot(jsub.exp, aes(y = rowcoord, x = colcoord, label = indx)) + geom_label()  + 
  theme_bw() + theme(aspect.ratio=2/3, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c() + 
  scale_y_continuous(trans = "reverse", breaks = seq(from = 0, to = 16, by = 4))  + 
  scale_x_continuous(breaks = seq(from = 0, to = 24, by = 4)) + 
  xlab("Column") + ylab("Row") + facet_wrap(~experi)


ggplot(jsub.sum, aes(y = rowcoord, x = colcoord, color = log2(totalcounts), size = log2(totalcounts))) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=2/3, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c() + 
  scale_y_continuous(trans = "reverse", breaks = seq(from = 0, to = 16, by = 4))  + 
  scale_x_continuous(breaks = seq(from = 0, to = 24, by = 4)) + 
  xlab("Column") + ylab("Row")  + facet_wrap(~experi)

ggplot(jsub.sum, aes(y = rowcoord, x = colcoord, color = log2(chromocounts), size = log2(chromocounts))) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=2/3, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c() + 
  scale_y_continuous(trans = "reverse", breaks = seq(from = 0, to = 16, by = 4))  + 
  scale_x_continuous(breaks = seq(from = 0, to = 24, by = 4)) + 
  xlab("Column") + ylab("Row") + facet_wrap(~experi)

ggplot(jsub.sum, aes(y = rowcoord, x = colcoord, color = log2(spikeincounts), size = log2(spikeincounts))) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=2/3, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c() + 
  scale_y_continuous(trans = "reverse", breaks = seq(from = 0, to = 16, by = 4))  + 
  scale_x_continuous(breaks = seq(from = 0, to = 24, by = 4)) + 
  xlab("Column") + ylab("Row") + facet_wrap(~experi)


# Summarize by column  ----------------------------------------------------

jsub.sum.bycol <- jsub.sum %>%
  group_by(experi, colcoord) %>%
  summarise(chromocounts = mean(chromocounts))

jsub.sum.bycol2 <- jsub.sum %>% ungroup() %>% 
  filter(chromocounts < 500000) %>%
  mutate(colcoord2 = factor(as.character(colcoord), levels = as.character(sort(unique(colcoord)))),
         colcoord.by4 = floor((colcoord - 1) / 4))

ggplot(jsub.sum.bycol, aes(y = log2(chromocounts), x = colcoord, color = log2(chromocounts), size = log2(chromocounts))) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=2/3, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c() + 
  scale_x_continuous(breaks = seq(from = 0, to = 24, by = 4)) + 
  facet_wrap(~experi) + ggtitle("Counts by column averages")

ggplot(jsub.sum.bycol2, aes(y = log2(chromocounts), x = colcoord2, fill = colcoord.by4)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=2/3, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c() + 
  facet_wrap(~experi) + ggtitle("Counts by column") + 
  xlab("Column") + 
  scale_fill_viridis_c()

ggplot(jsub.sum.bycol2, aes(y = chromocounts, x = colcoord2, fill = colcoord.by4)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=2/3, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c() + 
  facet_wrap(~experi) + ggtitle("Counts by column") + 
  xlab("Column") + 
  scale_fill_viridis_c()

ggplot(jsub.sum.bycol2 %>% filter(colcoord > 4), aes(y = log2(chromocounts), x = colcoord2)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=2/3, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c() + 
  facet_wrap(~experi) + ggtitle("Counts by column filt out emptys") + 
  xlab("Column")

ggplot(jsub.sum.bycol2 %>% filter(colcoord > 4), 
       aes(y = log2(chromocounts), x = colcoord2, fill = colcoord.by4)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=2/3, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c() + 
  facet_wrap(~experi) + ggtitle("Counts by column filt out emptys zoom") + 
  xlab("Column") + 
  scale_y_continuous(breaks = seq(10, 20, by = 1)) + 
  coord_cartesian(ylim = c(10, 20)) + 
  scale_fill_viridis_c()

ggplot(jsub.sum.bycol2 %>% filter(colcoord > 4), 
       aes(y = chromocounts, x = colcoord2, fill = colcoord.by4)) + geom_boxplot() + 
  theme_bw() + theme(aspect.ratio=2/3, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c() + 
  facet_wrap(~experi) + ggtitle("Counts by column filt out emptys zoom") + 
  xlab("Column") + 
  scale_fill_viridis_c()

# dev.off()


# Fraction counts on chromosome have less variabilllity in cells t --------

# show boxplots but by concentrations 

# create spikein? 

# conc2row <- 2^seq(from = log2(350), to = log2(50000), by = 1)
conc.vec <- rev(c(350, 700, 1500, 3000, 6000, 12000, 25000, 50000))
conc.long.vec <- rep(conc.vec, each = 2)
rows.vec <- seq(16)

ncells.vec <- seq(from = 0, to = 5)
ncells.long.vec <- rep(ncells.vec, each = 4)

cols.vec <- seq(24)

row2conc <- hash::hash(rows.vec, conc.long.vec)
col2ncells <- hash::hash(cols.vec, ncells.long.vec)

jsub.sum <- jsub.sum %>%
  rowwise() %>%
  mutate(spikeinconc = row2conc[[as.character(rowcoord)]],
         ncells = col2ncells[[as.character(colcoord)]])

jsub.sum$spikeinconcFactor <- factor(jsub.sum$spikeinconc, levels = unique(sort(jsub.sum$spikeinconc)))

# plot ratio of chromocount and spikein counts
ggplot(jsub.sum, aes(x = spikeinconcFactor, y = chromocounts / spikeincounts)) + 
  geom_boxplot() + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


mincounts <- 800
for (jprep in jconcs.lst){
  
  m0 <- ggplot(jsub.sum %>% filter(conc == jprep), aes(x = spikeinconcFactor, y = chromocounts / spikeincounts)) + 
    geom_boxplot() + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
    facet_grid(mark ~ ncells)  + ggtitle(paste("Preparation:", jprep))
  
  m1 <- ggplot(jsub.sum %>% filter(conc == jprep), aes(x = spikeinconcFactor, y = log2(chromocounts / spikeincounts))) + 
    geom_boxplot() + geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    facet_grid(mark ~ ncells) + ggtitle(paste("Preparation:", jprep))
  
  m1.filt <- ggplot(jsub.sum %>% filter(conc == jprep & chromocounts > mincounts), aes(x = spikeinconcFactor, y = log2(chromocounts / spikeincounts))) + 
    geom_boxplot() + geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    facet_grid(mark ~ ncells) + ggtitle(paste("Preparation:", jprep, "MinChromoCount >", mincounts))
  
  m2 <- ggplot(jsub.sum %>% filter(ncells > 0 & conc == jprep), aes(x = spikeinconcFactor, y = log2(chromocounts / spikeincounts))) + 
    geom_boxplot() + geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    facet_grid(mark ~ ncells) + ggtitle(paste("Preparation:", jprep))
  
  m2.filt <- ggplot(jsub.sum %>% filter(ncells > 0 & chromocounts > mincounts & conc == jprep), aes(x = spikeinconcFactor, y = log2(chromocounts / spikeincounts))) + 
    geom_boxplot() + geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
    facet_grid(mark ~ ncells) + ggtitle(paste("Preparation:", jprep, "MinChromoCount >", mincounts))
  
  
  print(m0)
  print(m1)
  print(m1.filt)
  print(m2)
  print(m2.filt)
}


ggplot(jsub.sum, aes(y = rowcoord, x = colcoord, color = ncells, size = ncells)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=2/3, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c() + 
  scale_y_continuous(trans = "reverse", breaks = seq(from = 0, to = 16, by = 4))  + 
  scale_x_continuous(breaks = seq(from = 0, to = 24, by = 4)) + 
  xlab("Column") + ylab("Row") + facet_wrap(~experi)

ggplot(jsub.sum, aes(y = rowcoord, x = colcoord, color = log2(spikeinconc), size = ncells)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=2/3, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  scale_color_viridis_c() + 
  scale_y_continuous(trans = "reverse", breaks = seq(from = 0, to = 16, by = 4))  + 
  scale_x_continuous(breaks = seq(from = 0, to = 24, by = 4)) + 
  xlab("Column") + ylab("Row") + facet_wrap(~experi)


# Calculate fold changes global -------------------------------------------

# fit linear model

print(unique(jsub.sum$spikeinconc))

input.dat <- subset(jsub.sum, spikeinconc == 50000 & conc == "37U" & mark == "H3K4me3" & ncells > 0)



spikeinconc.vec <- sort(unique(jsub.sum$spikeinconc))
conc.vec <- c("37U", "75U")
jmarks.vec <- c("H3K4me3", "H3K27me3")

for (jspikeinconc in spikeinconc.vec){
  for (jconc in conc.vec){
    for (jmark in jmarks.vec){
      input.dat <- subset(jsub.sum, spikeinconc == jspikeinconc & conc == jconc & mark == jmark & ncells > 0) %>%
        rowwise() %>%
        mutate(ncells.lin = ncells,
               ncells = log(ncells))  # makes fits linear
      # plot fits
      f1 <- FitNormCountsToNcells.lm(input.dat, return.fit.obj = TRUE)
      f2 <- FitNormCountsToNcells.glm(input.dat, return.fit.obj = TRUE)
      # f2.ci <- confint.default(f2)
      # jslope.glm.ln <- coefficients(f2)[["ncells"]]
      jpred1 <- data.frame(ypred = predict(f1, input.dat, se.fit = FALSE), ncells = input.dat$ncells, chromocounts = input.dat$chromocounts, spikeinconc = input.dat$spikeinconc, spikeincounts = input.dat$spikeincounts, stringsAsFactors = FALSE)
      jpred2 <- data.frame(ypred.unadj = predict(f2, input.dat, se.fit = FALSE), ncells = input.dat$ncells, chromocounts = input.dat$chromocounts, spikeinconc = input.dat$spikeinconc, spikeincounts = input.dat$spikeincounts, stringsAsFactors = FALSE) %>%
        mutate(ypred = ypred.unadj - log(spikeincounts))
      m.check1 <- ggplot(input.dat, 
             aes(x = ncells, y = log(chromocounts / spikeincounts))) + geom_point()  + 
        theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        geom_line(data = jpred1, mapping = aes(x = ncells, y = ypred)) + 
        ggtitle(paste("Lm fit: SpikeInMole:", jspikeinconc, jconc, jmark)) + 
        xlab("log(ncells)")
      print(m.check1)
      m.check2 <- ggplot(input.dat, 
             aes(x = ncells, y = log(chromocounts / spikeincounts))) + geom_point()  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        geom_line(data = jpred2, mapping = aes(x = ncells, y = ypred)) + 
        ggtitle(paste("GLM fit: SpikeInMole:", jspikeinconc, jconc, jmark)) + 
        xlab("log(ncells)")
      print(m.check2)
      # plot linear
      
      # jfit.l1 <- lm(formula = chromocounts ~ ncells.lin, data = input.dat)
      # jfit.l2 <- lm(formula = chromocounts / spikeincounts ~ ncells.lin, data = input.dat)
      m.l1 <- ggplot(input.dat, aes(x = ncells.lin, y = chromocounts)) + 
        geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        ggtitle(paste("Linear fit on linear scale SpikeInMole:", jspikeinconc, jconc, jmark)) + 
        geom_smooth(method = "lm", se = FALSE) 
      m.l2 <- ggplot(input.dat, aes(x = ncells.lin, y = chromocounts / spikeincounts)) + 
        geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
        ggtitle(paste("Linear fit on linear scale SpikeInMole:", jspikeinconc, jconc, jmark)) + 
        geom_smooth(method = "lm", se = FALSE) 
    }
  }
}


jsub.fits.lm <- jsub.sum %>% 
  filter(ncells > 0) %>%
  group_by(spikeinconc, conc, mark) %>%
  do(FitNormCountsToNcells.lm(.)) %>%
  ungroup() %>%
  mutate(spikeinconc.fct = factor(spikeinconc, levels = unique(sort(spikeinconc))))

jsub.fits.lm.naive <- jsub.sum %>% 
  filter(ncells > 0) %>%
  group_by(spikeinconc, conc, mark) %>%
  do(FitNormCountsToNcells.lm.naive(.)) %>%
  ungroup() %>%
  mutate(spikeinconc.fct = factor(spikeinconc, levels = unique(sort(spikeinconc))))

jsub.fits.glm <- jsub.sum %>% 
  filter(ncells > 0) %>%
  group_by(spikeinconc, conc, mark) %>%
  do(FitNormCountsToNcells.glm(.)) %>%
  ungroup() %>% 
  mutate(spikeinconc.fct = factor(spikeinconc, levels = unique(sort(spikeinconc))))


ggplot(jsub.fits.lm, aes(x = spikeinconc.fct, y = slope, ymin = slope - 2 * slope.se, ymax = slope + 2 * slope.se)) + geom_point() + 
  geom_errorbar() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_grid(mark ~ conc) + geom_hline(yintercept = 1, linetype = "dotted")  + 
  ggtitle("Fits on Total Counts: LM on chromocounts / spikeincounts")

ggplot(jsub.fits.lm.naive, aes(x = spikeinconc.fct, y = slope, ymin = slope - 2 * slope.se, ymax = slope + 2 * slope.se)) + geom_point() + 
  geom_errorbar() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_grid(mark ~ conc) + geom_hline(yintercept = 1, linetype = "dotted") + 
  ggtitle("Fits on Total Counts: LM on chromocounts only")

ggplot(jsub.fits.glm, aes(x = spikeinconc.fct, y = slope, ymin = slope - 2 * slope.se, ymax = slope + 2 * slope.se)) + geom_point() + 
  geom_errorbar() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_grid(mark ~ conc) + geom_hline(yintercept = 1, linetype = "dotted")  + 
  ggtitle("Fits on Total Counts: GLM on chromocounts / spikeincounts")


# Which spikein concentrations have most linear increase?  ----------------

# input.dat2 <- subset(jsub.sum, spikeinconc == 50000 & conc == "37U" & mark == "H3K4me3" & ncells > 0)
input.dat2 <- subset(jsub.sum, conc == "37U" & mark == "H3K4me3")
ggplot(input.dat2, aes(x = log2(spikeinconc), y = log2(spikeincounts))) + 
  geom_point(alpha = 0.25)  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~ncells) + geom_smooth(method = "lm", se = FALSE)



dat.rz.merge <- left_join(dat.rz, jsub.sum)

ggplot(dat.rz.merge, aes(x = chromocounts, y = TA.frac, color = log10(spikeinconc))) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() +
  scale_color_viridis_c() + 
  facet_grid(mark ~ conc) + 
  geom_vline(xintercept = mincounts)


dev.off()

# Save objects  -----------------------------------------------------------

saveRDS(jsub.sum, file = outrds)

# 
# 
# # More downstream achecks -------------------------------------------------
# 
# 
# m0 <- ggplot(jsub.sum %>% filter(conc == "37U"), aes(x = spikeinconcFactor, y = chromocounts / spikeincounts)) + 
#   geom_boxplot() + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
#   facet_grid(mark ~ ncells)  + ggtitle(paste("Preparation:", jprep))
# 
# m1 <- ggplot(jsub.sum %>% filter(conc == "37U"), aes(x = spikeinconcFactor, y = log2(chromocounts / spikeincounts))) + 
#   geom_boxplot() + geom_point() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
#   facet_grid(mark ~ ncells) + ggtitle(paste("Preparation:", jprep))
# 
# m2 <- ggplot(jsub.sum %>% filter(ncells > 0 & conc == "37U"), aes(x = spikeinconcFactor, y = log2(chromocounts / spikeincounts))) + 
#   geom_boxplot() + geom_point() + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
#   facet_grid(mark ~ ncells) + ggtitle(paste("Preparation:", jprep))
# 
# print(m0)
# print(m1)
# print(m2)


# what causes variability?  -----------------------------------------------

jprep <- "75U"
m1 <- ggplot(jsub.sum %>% filter(conc == jprep), aes(x = spikeinconcFactor, y = log2(chromocounts / spikeincounts))) + 
  geom_boxplot() + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_grid(mark ~ ncells) + ggtitle(paste("Preparation:", jprep))

m1.top <- ggplot(jsub.sum %>% filter(conc == jprep), aes(x = spikeinconcFactor, y = log2(chromocounts))) + 
  geom_boxplot() + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_grid(mark ~ ncells) + ggtitle(paste("Preparation:", jprep))

m1.bottom <- ggplot(jsub.sum %>% filter(conc == jprep), aes(x = spikeinconcFactor, y = log2(spikeincounts))) + 
  geom_boxplot() + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_grid(mark ~ ncells) + ggtitle(paste("Preparation:", jprep))

print(m1)

multiplot(m1.top, m1.bottom, cols = 1)


# Those bad chromocounts are just bad samples? ----------------------------

dat.rz.merge <- left_join(dat.rz, jsub.sum)

ggplot(dat.rz.merge, aes(x = chromocounts, y = TA.frac, color = log10(spikeinconc))) + 
  geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_x_log10() +
  scale_color_viridis_c() + 
  facet_grid(mark ~ conc) + 
  geom_vline(xintercept = mincounts)


# x <- rep(seq(5), each = 3)
# y <- 100 * x 
# y.noisy <- y + rnorm(n = length(y), mean = 0, sd = 2)
# 
# plot(x, y.noisy)
# plot(log(x), log(y.noisy))
# 
# jdat <- data.frame(x = x, y.noisy = y.noisy, stringsAsFactors = FALSE)
# lm(formula = y.noisy ~ x, data = jdat)
# lm(formula = log(y.noisy) ~ log(x), data = jdat)
# lm(formula = log2(y.noisy) ~ log2(x), data = jdat)
# lm(formula = log(y.noisy) ~ x, data = jdat)
  

