# Jake Yeung
# Date of Creation: 2020-08-10
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/predict_K562_from_spikeins.R
# description

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)


# Load different spikeins measure counts in each chromosome ------------------------------------------------

outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins/genomewide_analysis_spikein_vs_naive.", Sys.Date(), ".pdf")

# pdf(pdfout, useDingbats = FALSE)
pdf(outpdf, width = 1020/72, height = 815/72, useDingbats = FALSE)

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



# Predict number of cells from total counts?  -----------------------------

# remove bad cells


# ggplot(dat.filt.long %>% filter(chromocounts >= 1000), aes(x = log2(chromocounts / spikeincounts))) + 
#   geom_density() + facet_grid(conc~spikeinconcFactor) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


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



jprep <- "75U"
jprep <- "37U"
mincounts <- 1000
m1 <- ggplot(jsub.sum %>% filter(conc == jprep), aes(x = spikeinconcFactor, y = log2(chromocounts / spikeincounts))) + 
  geom_boxplot() + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
  facet_grid(mark ~ ncells) + ggtitle(paste("Preparation:", jprep))
print(m1)


ggplot(jsub.sum %>% filter(conc == jprep), 
       aes(x = log2(chromocounts / spikeincounts), fill = as.character(ncells))) + 
  geom_density(alpha = 0.25) + 
  scale_fill_viridis_d() + 
  facet_grid(mark ~ spikeinconcFactor) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(jsub.sum %>% filter(conc == jprep & chromocounts > mincounts & ncells > 0), 
       aes(x = log2(chromocounts / spikeincounts), fill = as.character(ncells))) + 
  geom_density(alpha = 0.25) + 
  scale_fill_viridis_d() + 
  facet_grid(mark ~ spikeinconcFactor) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Preparation:", jprep)


ggplot(jsub.sum %>% filter(conc == jprep & chromocounts > mincounts & ncells > 0), 
       aes(x = log2(chromocounts / spikeincounts))) + 
  geom_density(alpha = 0.25) + 
  scale_fill_viridis_d() + 
  facet_grid(mark ~ spikeinconcFactor) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(jsub.sum %>% filter(conc == jprep & chromocounts > mincounts & ncells > 0), 
       aes(x = log2(chromocounts), fill = as.character(ncells))) + 
  geom_density(alpha = 0.25) + 
  scale_fill_viridis_d() + 
  facet_grid(mark ~ spikeinconcFactor) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Plot linear fit  --------------------------------------------------------

# jprep <- "37U"
# jprep <- "75U"

jpreps <- c("37U", "75U")

for (jprep in jpreps){
  
  m1 <- ggplot(jsub.sum %>% filter(conc == jprep & chromocounts > mincounts & ncells > 0), 
         aes(y = log2(chromocounts / spikeincounts), x = log(ncells))) +
    geom_point() + 
    scale_fill_viridis_d() + 
    facet_grid(mark ~ spikeinconcFactor) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(paste0("Preparation:", jprep)) + 
    geom_smooth(method = "lm", se = TRUE)
  
  
  m2 <- ggplot(jsub.sum %>% filter(conc == jprep & chromocounts > mincounts & ncells > 0), 
         aes(y = log2(chromocounts), x = log(ncells))) +
    geom_point() + 
    scale_fill_viridis_d() + 
    facet_grid(mark ~ spikeinconcFactor) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(paste0("Preparation:", jprep)) + 
    geom_smooth(method = "lm", se = TRUE)
  
  print(m1)
  print(m2)
  
}


# do the fits

jsub.sum.filt <- jsub.sum %>% filter(chromocounts > mincounts & ncells > 0)

# # remove outlier
# # jsub.sum.filt <- subset(jsub.sum.filt, !(spikeinconc == 6000 & ncells == 2 & mark == "H3K4me3"))
# jsubsub <- jsub.sum.filt %>% filter(spikeinconc == 350 & mark == "H3K27me3" & conc == jprep) 
# 
# input.dat <- jsubsub %>%
#   rowwise() %>%
#   mutate(genecounts = chromocounts,
#          ncells = log(ncells))
# 
# input.raw.dat <- jsubsub %>%
#   rowwise() %>%
#   mutate(genecounts = chromocounts,
#          chromocounts = 1, 
#          ncells = log(ncells))
# 
# ggplot(input.dat, aes(x = ncells, y = log2(genecounts  / spikeincounts))) +
#   geom_point() + 
#   xlab("log(ncells)") + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# ggplot(input.raw.dat, aes(x = ncells, y = log2(genecounts))) +
#   geom_point() + 
#   xlab("log(ncells)") + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#  


# Plot for all  -----------------------------------------------------------

jfits.all.spike <- jsub.sum.filt %>%
  group_by(spikeinconcFactor, mark, conc) %>%
  mutate(genecounts = chromocounts,
         ncells = log(ncells)) %>%
  do(FitLmRowSpikeins(input.dat = ., return.fit.obj = FALSE, pseudocount = 0)) %>%
  mutate(type = "spikein")

jfits.all.naive <- jsub.sum.filt %>%
  group_by(spikeinconcFactor, mark, conc) %>%
  mutate(genecounts = chromocounts,
         chromocounts = 1,
         ncells = log(ncells)) %>%
  do(FitLmRowChromocounts(input.dat = ., return.fit.obj = FALSE, pseudocount = 0)) %>%
  mutate(type = "naive")
  
ggplot(rbind(jfits.all.spike, jfits.all.naive), aes(x = slope.ln, fill = type)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  ggtitle("Slopes result from each spikeInConcentration, preparation, and mark") + 
  xlab("Slope")

ggplot(rbind(jfits.all.spike, jfits.all.naive), aes(x = slope.se.ln, fill = type)) + 
  geom_density(alpha = 0.25) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
  ggtitle("Slopes result from each spikeInConcentration, preparation, and mark") + 
  xlab("Standard Error of the Slope")

dev.off()


# Pretend it's two plates -------------------------------------------------


m1 <- ggplot(jsub.sum %>% filter(chromocounts > mincounts & ncells > 0), 
       aes(y = log2(chromocounts / spikeincounts), x = log(ncells), color = conc)) +
  geom_point(alpha = 0.25) + 
  scale_fill_viridis_d() + 
  facet_grid(mark ~ spikeinconcFactor) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_smooth(method = "lm", se = TRUE)
print(m1)  

m1 <- ggplot(jsub.sum %>% filter(chromocounts > mincounts & ncells > 0), 
       aes(y = log2(chromocounts), x = log(ncells), color = conc)) +
  geom_point(alpha = 0.25) + 
  facet_grid(mark ~ spikeinconcFactor) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_smooth(method = "lm", se = TRUE)
print(m1)  
