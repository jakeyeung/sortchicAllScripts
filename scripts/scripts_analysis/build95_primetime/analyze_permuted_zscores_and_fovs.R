# Jake Yeung
# Date of Creation: 2019-05-01
# File: ~/projects/scchic/scripts/scripts_analysis/build95_primetime/analyze_permuted_zscores_and_fovs.R
# Do analysis on FOV and Zscores permutations 

rm(list=ls())

library(dplyr)
library(ggplot2)
library(data.table)
library(ggrepel)

source("scripts/Rfunctions/BackgroundPermutationScripts.R")
source("scripts/Rfunctions/MaraDownstream.R")
source("scripts/Rfunctions/MatchCellNameToSample.R")

jmark <- "H3K4me1"

pdf("~/data/scchic/pdfs/background_model_motifs_MARA.pdf", useDingbats = FALSE)

# Load MARA  --------------------------------------------------------------

experis <- read.table("data/outputs/experinames.txt", stringsAsFactors = FALSE)
experihash <- hash(sapply(experis$V1, function(x) strsplit(x, "_")[[1]][[1]]), sapply(experis$V1, function(x) paste(strsplit(x, "_")[[1]][2:3], collapse="_")))
switch.rep.hash <- GetRepSwitchHash(experihash)
# add m1_S1 which is a bug from m1_S10
# fix a bug with GetTechRep which causes S1 to show up 
switch.rep.hash[["BM_H3K4me3_m1_S1"]] <- "BM_H3K4me3_m1_rep1"


indir.mara1 <- "/Users/yeung/data/scchic/from_cluster/mara_analysis_cluster_build95_CorrPeakFilt.withchr.cells_from_bin_analysis/H3K4me1/mara_output/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_H3K4me1.filt_0.99.center_TRUE_K50-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--K50/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_H3K4me1.filt_0.99.center_TRUE_K50"
# indir.mara2 <- "/Users/yeung/data/scchic/from_cluster/mara_analysis_cluster_build95_CorrPeakFilt.withchr.cells_from_bin_analysis/H3K4me3/mara_output/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_H3K4me3.filt_0.99.center_TRUE_K50-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--K50/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_H3K4me3.filt_0.99.center_TRUE_K50"

assertthat::assert_that(dir.exists(indir.mara1))
# assertthat::assert_that(dir.exists(indir.mara2))

mara.out <- LoadMARA(indir.mara1, rep.prefix = "", swap.tech.rep = switch.rep.hash, filesuffix = "")
# mara.out2 <- LoadMARA(indir.mara2, rep.prefix = "", swap.tech.rep = switch.rep.hash, filesuffix = "")

mara.outs <- list("H3K4me1" = mara.out)

# Load data ---------------------------------------------------------------

inf.fov <- "/Users/yeung/data/scchic/from_cluster/zscore_permute_summary_CorrPeakFilt/H3K4me1/FOV_permute_summary.txt.gz"

dat.fov <- fread(inf.fov, stringsAsFactors = FALSE, col.names = c("fov", "seed", "mark"))

marks <- c("H3K4me1")

inf.fov.real <- "/Users/yeung/data/scchic/from_cluster/mara_analysis_cluster_build95_CorrPeakFilt.withchr.cells_from_bin_analysis/H3K4me1/mara_output/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_H3K4me1.filt_0.99.center_TRUE_K50-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.bugfix.filt.0--K50/hiddenDomains_cellmin_0-cellmax_9999999-binarize_TRUE-BM_H3K4me1.filt_0.99.center_TRUE_K50/FOV"
fovs.real <- unlist(fread(inf.fov.real), use.names = FALSE)
fovs.real <- data.frame(fov.real = unlist(fovs.real), mark = marks[[1]])

dat.fov <- left_join(dat.fov, fovs.real)

pvals.long <- dat.fov %>%
  group_by(mark) %>%
  do(GetPvalFOV(., fovs.real = NULL, jprob = 0.9, show.plot = FALSE, return.pval.only = TRUE))

pvals.fov <- GetPvalFOV(dat.fov, fovs.real = NULL, jprob = 0.9, show.plot = FALSE, return.pval.only = FALSE)

m1 <- ggplot(pvals.fov$real.dat, aes(x = fov)) + geom_histogram(bins = 40) + theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle("FOV") + 
  xlab("FOV from Randomized Data")
m2 <- ggplot(pvals.fov$real.dat, aes(x = fov, y = log10.frac.more.than)) + geom_point() +
  xlab("FOV") + ylab("log10 Fraction Randomized Data >= FOV") + ggtitle("FOV") + 
  theme_bw(24) +
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = pvals.fov$fov.real, linetype = "dashed") + 
  expand_limits(y = ceiling(pvals.fov$log10pval)) + 
  geom_line(mapping = aes(x = fov, y = log10.frac.more.than), data = pvals.fov$pred.dat)
m3 <- ggplot(pvals.fov$real.dat, aes(x = fov, y = log10.frac.more.than)) + geom_point() +
  xlab("FOV") + ylab("log10 Fraction Randomized Data >= FOV") + ggtitle("FOV") + 
  theme_bw(24) + 
  theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_line(mapping = aes(x = fov, y = log10.frac.more.than), data = subset(pvals.fov$pred.dat, fov < max(pvals.fov$real.dat$fov)))
print(m1)
print(m2)
print(m3)


m1.fov <- ggplot(dat.fov, aes(x = fov)) + geom_histogram(bins = 40) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_vline(xintercept = fovs.real$fov.real[[1]])

m1.fov2 <- ggplot(dat.fov, aes(x = fov)) + geom_histogram(bins = 40) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Load zscore mats --------------------------------------------------------

jmarks <- c("H3K4me1")
names(jmarks) <- jmarks

inf.zscores <- lapply(jmarks, function(jmark) paste0("/Users/yeung/data/scchic/from_cluster/zscore_permute_summary_CorrPeakFilt/", jmark, "/zscore_permute_summary.txt.gz"))

dats <- lapply(inf.zscores, function(inf.zscore) fread(inf.zscore, header = FALSE, stringsAsFactors = FALSE))

zscore.long <- lapply(dats, function(dat){
  colnames(dat) <- c("motif", "zscore", "seed", "mark")
  return(dat)
})

zscore.long <- do.call(rbind, zscore.long)

print(head(zscore.long))
print(unique(zscore.long$motif))

zscore.long <- subset(zscore.long, mark != "exiting")
subset(zscore.long, grepl("Zscores", motif))

zscore.long$zscore <- as.numeric(zscore.long$zscore)

zscore.long.real <- lapply(jmarks, function(jmark){
  return(mara.outs[[jmark]]$zscores %>% mutate(mark = jmark))
}) %>% 
  bind_rows() %>%
  dplyr::rename(zscore.real = zscore)

zscore.long <- left_join(zscore.long, zscore.long.real)

zscore.long.real <- zscore.long.real %>%
  group_by(mark) %>%
  mutate(indx = seq(length(motif)),
         motif.lab = ifelse(zscore.real > 2, motif, NA))
  


# Plot CEBPB zscore versus background model  ------------------------------

# out <- GetPvalZscore(zscore.long %>% filter(mark == jmark & motif == jmotif), subset(mara.outs$H3K4me1$zscore, motif == jmotif)$zscore)

pvals.long <- zscore.long %>%
  group_by(mark, motif) %>%
  do(GetPvalZscore(., zscore.real = NULL, jprob = 0.9, show.plot = FALSE, return.pval.only = TRUE))

pvals.long <- pvals.long %>%
  group_by(mark) %>%
  arrange(log10pval) %>%
  mutate(indx = seq(length(motif))) %>%
  mutate(motif.lab = ifelse(log10pval < -35, motif, NA))

# plot top hits by mark
m.top <- ggplot(pvals.long, aes(x = log10pval)) + geom_histogram() + facet_wrap(~mark, ncol = 1) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# plot names of top hits?
m.top.name <- ggplot(pvals.long, aes(y = -log10pval, x = indx, label = motif.lab)) + geom_point() + geom_text_repel() + facet_wrap(~mark, ncol = 1) + 
  theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# plot real zscores
m.top.real.zscore <- ggplot(zscore.long.real, aes(x = indx, y = zscore.real, label = motif.lab)) + 
  geom_point() + geom_text_repel() + theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  xlab("Index") + ylab("Activity Zscore")

print(m.top)
print(m.top.name)
print(m.top.real.zscore)
 
# show individual hits
jmotifs <- c("Cebpb", "Ebf1", "Tal1", "Ets1", "Cebpa", "Foxc1")
for (jmotif in jmotifs){
  print(jmotif)
  pval <- GetPvalZscore(zscore.long %>% filter(motif == jmotif & mark == jmark), zscore.real = NULL, jprob = 0.9, show.plot = FALSE)
  m1 <- ggplot(pval$real.dat, aes(x = zscore)) + geom_histogram(bins = 40) + theme_bw(24) + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(jmotif) + 
    xlab("Zscore from Randomized Data")
  m2 <- ggplot(pval$real.dat, aes(x = zscore, y = log10.frac.more.than)) + geom_point() +
    xlab("Zscore") + ylab("log10 Fraction Randomized Data >= Zscore") + ggtitle(jmotif) + 
    theme_bw(24) +
    theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_vline(xintercept = pval$zscore.real, linetype = "dashed") + 
    expand_limits(y = ceiling(pval$log10pval)) + 
    geom_line(mapping = aes(x = zscore, y = log10.frac.more.than), data = pval$pred.dat)
  m3 <- ggplot(pval$real.dat, aes(x = zscore, y = log10.frac.more.than)) + geom_point() +
    xlab("Zscore") + ylab("log10 Fraction Randomized Data >= Zscore") + ggtitle(jmotif) + 
    theme_bw(24) + 
    theme(aspect.ratio=0.5, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_line(mapping = aes(x = zscore, y = log10.frac.more.than), data = subset(pval$pred.dat, zscore < max(pval$real.dat$zscore)))
  print(m1)
  print(m2)
  print(m3)
}


# show top hits
pvals.top <- pvals.long %>%
  group_by(mark) %>%
  dplyr::top_n(n = 10, wt = -log10pval) %>%
  arrange(log10pval)
print(dim(pvals.top))
print(split(pvals.top, f = pvals.top$mark))



dev.off()