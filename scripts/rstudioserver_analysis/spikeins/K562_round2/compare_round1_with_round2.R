# Jake Yeung
# Date of Creation: 2020-09-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/K562_round2/compare_round1_with_round2.R
# Compare Round1 vs Round2 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# Load spikeins -----------------------------------------------------------


hubprefix <- "/home/jyeung/hub_oudenaarden"
indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN4969/K562/tagged_bams/countTablesAndRZr1only_ByChromo.NewFilters")
infs.chromo <- list.files(indir, pattern = "K562-EtOH-.*.csv", full.names = TRUE)

jspikeinchromo <- "J02459.1"
jchromos <- paste("", seq(19), sep = "")

dat.chromos <- lapply(infs.chromo, function(inf){
  dat.filt.long <- GetChromoCounts(inf, spikeinchromo = jspikeinchromo, chromos.keep = jchromos)
}) %>%
  bind_rows() %>%
  filter(chromo == "1")

dat.spikeins.mat1 <- as.data.frame(subset(dat.chromos, select = c(samp, spikeincounts, chromocounts)))
rownames(dat.spikeins.mat1) <- dat.spikeins.mat1$samp

dat.spikeins.mat1$round <- "round1"

# 
# # Load spikeins round1 ----------------------------------------------------
# 
# inf.r1 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562/spikein_counts/spikein_counts_all.RData"
# load(inf.r1, v=T)
# 
# dat.spikeins.mat1 <- dat.spikeins.mat %>%
#   mutate(round = "round1")

# Load spikeins round 2 ---------------------------------------------------


# load data  ---------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

indir.chromo.g1 <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2/tagged_bams/G1_sorted/merged_bams/countTablesAndRZr1only_ByChromo.NewFilters")
indir.chromo.cc <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/VAN5039_K562_round2/tagged_bams/cellcycle_sorted/merged_bams/countTablesAndRZr1only_ByChromo.NewFilters")

assertthat::assert_that(dir.exists(indir.chromo.g1))
assertthat::assert_that(dir.exists(indir.chromo.cc))

infs.chromo.g1 <- list.files(indir.chromo.g1, pattern = "K562-EtOH-.*.csv", full.names = TRUE)
infs.chromo.cc <- list.files(indir.chromo.cc, pattern = "K562-EtOH-.*.csv", full.names = TRUE)


jspikeinchromo <- "J02459.1"
jchromos <- paste("", seq(19), sep = "")

dat.chromos <- lapply(c(infs.chromo.g1, infs.chromo.cc), function(inf){
  dat.filt.long <- GetChromoCounts(inf, spikeinchromo = jspikeinchromo, chromos.keep = jchromos)
}) %>%
  bind_rows()

# dat.chromos.g1 <- lapply(infs.chromo.g1, function(inf){
#   dat.filt.long <- GetChromoCounts(inf, spikeinchromo = jspikeinchromo, chromos.keep = jchromos)
# }) %>%
#   bind_rows()
# 
# dat.chromos.cc <- lapply(infs.chromo.cc, function(inf){
#   dat.filt.long <- GetChromoCounts(inf, spikeinchromo = jspikeinchromo, chromos.keep = jchromos)
# }) %>%
#   bind_rows()

chromocounts <- subset(dat.chromos, chromo == "1", select = c(samp, chromocounts, spikeincounts))

dat.spikeins.mat2 <- as.data.frame(chromocounts)

# add spikeins
rownames(dat.spikeins.mat2) <- dat.spikeins.mat2$samp
dat.spikeins.mat2$round <- "round2"


# 
# inf.r2 <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/quality_control_spikeins_K562.round2/dat_spikeins_all.round2.RData"
# load(inf.r2, v=T)
# 
# dat.spikeins.mat2 <- dat.spikeins.mat %>%
#   mutate(round = "round2")

cells.common <- intersect(dat.spikeins.mat1$samp, dat.spikeins.mat2$samp)

# Compare spikeincounts ---------------------------------------------------

dat.spikeins.mat.merged <- rbind(dat.spikeins.mat1, dat.spikeins.mat2) %>%
  filter(samp %in% cells.common) %>%
  rowwise() %>%
  mutate(log2ratio = log2(chromocounts / spikeincounts),
         experi = ClipLast(samp, jsep = "_"))


spikeins.wide <- dcast(data = dat.spikeins.mat.merged, formula = samp + experi ~ round, value.var = "spikeincounts")
chromos.wide <- dcast(data = dat.spikeins.mat.merged, formula = samp + experi ~ round, value.var = "chromocounts")
log2ratio.wide <- dcast(data = dat.spikeins.mat.merged, formula = samp + experi~ round, value.var = "log2ratio")

log2ratio.wide$mark <- sapply(log2ratio.wide$experi, function(x) strsplit(x, "-")[[1]][[3]])
log2ratio.wide$plate <- sapply(log2ratio.wide$experi, function(x){
  xtmp <- strsplit(x, "-")[[1]]
  xtmp <- xtmp[4:length(xtmp)]
  xout <- paste(xtmp, collapse = "-")
})

ggplot(spikeins.wide, aes(x = round1, y = round2, color = experi)) +  
  geom_point()  + 
  ggtitle("Number of spikein cuts") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(spikeins.wide, aes(x = round1, y = round2, color = plate)) +  
  geom_point()  + 
  ggtitle("Number of spikein cuts") + 
  facet_wrap(~mark) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(chromos.wide, aes(x = round1, y = round2, color = experi)) + 
  geom_point() + 
  ggtitle("Number of chromo cuts") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(log2ratio.wide, aes(x = round1, y = round2, color = experi)) + 
  geom_point(alpha = 0.25) + 
  ggtitle("log2(chromo / spike) cuts") + 
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  scale_x_continuous(breaks = seq(0, 8)) + 
  scale_y_continuous(breaks = seq(0, 8)) 

