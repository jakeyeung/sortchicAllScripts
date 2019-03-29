# Jake Yeung
# Date of Creation: 2019-03-24
# File: ~/projects/scchic/scripts/scripts_analysis/pseudobulk_analysis/plot_correlations.R
# Plot correlations

library(ggplot2)
library(dplyr)
library(reshape2)



# Functions ---------------------------------------------------------------




ComparePublicLog <- function(inf, thres = 0.995, lab = "MyLabel"){
  dat <- data.table::fread(inf, sep = "\t")
  
  datref <- dat[, 1]
  colnames(datref) <- c("chip")
  datcompare <- dat[, 2:ncol(dat)]
  
  datref$peakid <- seq(nrow(datref))
  datcompare$peakid <- seq(nrow(datcompare))
  
  # plot 1 vs all
  datcompare.long <- melt(datcompare, variable.name = "Sample", value.name = "chic", id.vars = "peakid")
  
  dat.long <- dplyr::left_join(datcompare.long, datref)
  
  dat.long.filt <- dat.long %>%
    filter(!is.nan(chip)) %>%
    group_by(Sample) %>%
    filter(chic < quantile(chic, probs = thres, na.rm = TRUE)) %>%
    mutate(compare = lab)
  
  
  dat.cors <- dat.long.filt %>%
    group_by(Sample) %>%
    summarise(corr.pears = cor(chip, chic, method = "pearson"),
              corr.spear = cor(chip, chic, method = "spearman")) %>%
    mutate(compare = lab)
  # remoev outliers?
  
  # m <- ggplot(dat, aes(x = H3K4me1_MatBcell_rep1_reltoinput.bw, y = H3K4me1_cluster_1.log1p.bw)) + geom_point() + 
  #   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m <- ggplot(dat.long.filt, aes(x = chip, y = chic)) + 
    geom_hex(bins = 25) +
    # geom_point(alpha = 0.01) + 
    # stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~Sample)
  return(list(dat.long.filt = dat.long.filt, dat.cors = dat.cors, m = m))
}

ComparePublicLinear <- function(inf, thres = 0.995, lab = "MyLabel"){
  dat <- data.table::fread(inf, sep = "\t")
  
  datref <- dat[, 1]
  colnames(datref) <- c("chip")
  datcompare <- dat[, 2:ncol(dat)]
  
  datref$peakid <- seq(nrow(datref))
  datcompare$peakid <- seq(nrow(datcompare))
  
  # plot 1 vs all
  datcompare.long <- melt(datcompare, variable.name = "Sample", value.name = "chic", id.vars = "peakid")
  
  dat.long <- dplyr::left_join(datcompare.long, datref)
  
  # remove NaNs
  
  dat.long.filt <- dat.long %>%
    filter(!is.nan(chip)) %>%
    group_by(Sample) %>%
    filter(log(chic) < quantile(log(chic), probs = thres, na.rm = TRUE) & log(chic) > quantile(log(chic), probs = 1 - thres)) %>%
    filter(log(chip) < quantile(log(chip), probs = thres, na.rm = TRUE) & log(chip) > quantile(log(chip), probs = 1 - thres)) %>%
    mutate(compare = lab)
  
  # dat.long.filt <- subset(dat.long.filt, !is.nan(chip))
  
  dat.cors <- dat.long.filt %>%
    group_by(Sample) %>%
    summarise(corr.pears = cor(chip, chic, method = "pearson"),
              corr.spear = cor(chip, chic, method = "spearman"),
              corr.pears.log = cor(log(chip), log(chic), method = "pearson"),
              corr.spear.log = cor(log(chip), log(chic), method = "spearman")) %>%
    mutate(compare = lab)
  # remoev outliers?
  
  # m <- ggplot(dat, aes(x = H3K4me1_MatBcell_rep1_reltoinput.bw, y = H3K4me1_cluster_1.log1p.bw)) + geom_point() + 
  #   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m <- ggplot(dat.long.filt, aes(x = log10(chip), y = log10(chic))) + 
    geom_hex(bins = 25) +
    # geom_point(alpha = 0.01) + 
    # stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~Sample)
  return(list(dat.long.filt = dat.long.filt, dat.cors = dat.cors, m = m))
}

CompareAcrossMarks <- function(inf2, jmark, clstr, thres = 0.995, lab = "MyLabel"){
  dat.multi <- data.table::fread(inf2, sep = "\t")
  # compare erythryroblasts
  datref.multi <- dat.multi[, ..clstr]
  colnames(datref.multi) <- c("chic_ref")
  cols.compare <- grepl(jmark, colnames(dat.multi))
  datcompare.multi <- dat.multi[, ..cols.compare]
  
  datref.multi$peakid <- seq(nrow(datref.multi))
  datcompare.multi$peakid <- seq(nrow(datcompare.multi))
  
  datcompare.multi.long <- melt(datcompare.multi, id.vars = "peakid", variable.name = "Sample", value.name = "chic")
  
  dat.multi.long <- dplyr::left_join(datcompare.multi.long, datref.multi)
  
  # remove outliers
  dat.multi.long.filt <- dat.multi.long %>%
    group_by(Sample) %>%
    filter(chic < quantile(chic, probs = thres) & chic > quantile(chic, probs = 1 - thres)) %>%
    filter(chic_ref < quantile(chic_ref, probs = thres) & chic_ref > quantile(chic_ref, probs = 1 - thres)) %>%
    mutate(compare = lab)
  
  m <- ggplot(dat.multi.long.filt, aes(x = chic_ref, y = chic)) + geom_point(alpha = 0.25) + 
    scale_x_log10() + scale_y_log10() + facet_wrap(~Sample) + theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  dat.cors <- dat.multi.long.filt %>%
    group_by(Sample) %>%
    summarise(corr.pears = cor(chic, chic_ref, method = "pearson"),
              corr.spear = cor(chic, chic_ref, method = "spearman"),
              corr.pears.log2 = cor(log2(chic), log2(chic_ref), method = "pearson"),
              corr.spear.log2 = cor(log2(chic), log2(chic_ref), method = "spearman")) %>%
    
    mutate(compare = lab)
  # print(dat.cors)
  return(list(dat.multi.long.filt = dat.multi.long.filt, dat.cors = dat.cors, m = m))
}

# With log ----------------------------------------------------------------


jmark <- "H3K4me3"

jmark <-  "H3K4me1"

jmark <- "H3K27me3"
jmark <- "H3K27me3"
inf <- "/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/H3K4me1_Neu_comparison_InputNorm.tsv"
inf <- "/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/H3K4me3_Erythrobl_rep1_comparison_InputNorm.tsv"
inf <- "/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/H3K4me3_Megakaryo_rep1_comparison.tsv"
inf <- "/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/H3K4me1_Megakaryo_rep1_comparison.tsv"
inf <- "/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/H3K4me1_MatBcell_rep1_comparison_InputNorm.tsv"
# inf <- "/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/H3K4me1_MatBcell_rep2_comparison_InputNorm.tsv"

# inf <- "/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/H3K4me1_Neu_comparison_InputNorm.tsv"

inf <- paste0("/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/", jmark, "_MatBcell_rep1_comparison_InputNorm.tsv")
inf <- paste0("/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/", jmark, "_Neu_comparison_InputNorm.tsv")


assertthat::assert_that(file.exists(inf))



dat <- data.table::fread(inf, sep = "\t")

datref <- dat[, 1]
colnames(datref) <- c("chip")
datcompare <- dat[, 2:ncol(dat)]

datref$peakid <- seq(nrow(datref))
datcompare$peakid <- seq(nrow(datcompare))

# plot 1 vs all
datcompare.long <- melt(datcompare, variable.name = "Sample", value.name = "chic", id.vars = "peakid")

dat.long <- dplyr::left_join(datcompare.long, datref)

dat.long.filt <- dat.long %>%
  group_by(Sample) %>%
  filter(chic < quantile(chic, probs = thres, na.rm = TRUE))

dat.long.filt <- subset(dat.long.filt, !is.nan(chip))

dat.cors <- dat.long.filt %>%
  group_by(Sample) %>%
  summarise(corr.pears = cor(chip, chic, method = "pearson"),
            corr.spear = cor(chip, chic, method = "spearman")) %>%
  mutate(compare = lab)
# remoev outliers?

# m <- ggplot(dat, aes(x = H3K4me1_MatBcell_rep1_reltoinput.bw, y = H3K4me1_cluster_1.log1p.bw)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m <- ggplot(dat.long.filt, aes(x = chip, y = chic)) + 
  geom_hex(bins = 25) +
  # geom_point(alpha = 0.01) + 
  # stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~Sample)



# Without log -------------------------------------------------------------

jmark <- "H3K4me3"
jmark <- "H3K4me1"

jmark <- "H3K27me3"
thres <- 0.95
inf <- paste0("/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/", jmark, "_MatBcell_MatBcell_comparison.tsv")
# inf <- paste0("/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/", jmark, "_Neu_comparison.tsv")
assertthat::assert_that(file.exists(inf))

dat <- data.table::fread(inf, sep = "\t")

datref <- dat[, 1]
colnames(datref) <- c("chip")
datcompare <- dat[, 2:ncol(dat)]

datref$peakid <- seq(nrow(datref))
datcompare$peakid <- seq(nrow(datcompare))

# plot 1 vs all
datcompare.long <- melt(datcompare, variable.name = "Sample", value.name = "chic", id.vars = "peakid")

dat.long <- dplyr::left_join(datcompare.long, datref)

# remove NaNs

dat.long.filt <- dat.long %>%
  filter(!is.nan(chip)) %>%
  group_by(Sample) %>%
  filter(log(chic) < quantile(log(chic), probs = thres, na.rm = TRUE) & log(chic) > quantile(log(chic), probs = 1 - thres)) %>%
  filter(log(chip) < quantile(log(chip), probs = thres, na.rm = TRUE) & log(chip) > quantile(log(chip), probs = 1 - thres))

# dat.long.filt <- subset(dat.long.filt, !is.nan(chip))

dat.cors <- dat.long.filt %>%
  group_by(Sample) %>%
  summarise(corr.pears = cor(chip, chic, method = "pearson"),
            corr.spear = cor(chip, chic, method = "spearman"))
# remoev outliers?

print(dat.cors)

# m <- ggplot(dat, aes(x = H3K4me1_MatBcell_rep1_reltoinput.bw, y = H3K4me1_cluster_1.log1p.bw)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

m <- ggplot(dat.long.filt, aes(x = log10(chip), y = log10(chic))) + 
  geom_hex(bins = 25) +
  # geom_point(alpha = 0.01) + 
  # stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~Sample)
print(m)



# Correlate mark between mark ---------------------------------------------

jmark <- "H3K4me3"
jmark <- "H3K27me3"

jmark <- "H3K9me3"
# inf2 <- "/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/H3K9me3_Multimark_H3K4me1_vs_H3K9me3_comparison.tsv"
inf2 <- paste0("/Users/yeung/data/scchic/public_data/bigwig_compares_output_tables/", jmark, "_Multimark_H3K4me1_vs_", jmark, "_comparison.tsv")

assertthat::assert_that(inf2)

is.log <- FALSE

dat.multi <- data.table::fread(inf2, sep = "\t")

# compare erythryroblasts
clstr <- 1  # Bcells

clstr <- 2  # Bcells upstream

clstr <- 3  # Neutrophils

clstr <- 9  # Neutrohils upstream 

clstr <- 4  # NK cells??

datref.multi <- dat.multi[, ..clstr]
print(head(datref.multi))
colnames(datref.multi) <- c("chic_ref")
cols.compare <- grepl(jmark, colnames(dat.multi))
datcompare.multi <- dat.multi[, ..cols.compare]

datref.multi$peakid <- seq(nrow(datref.multi))
datcompare.multi$peakid <- seq(nrow(datcompare.multi))

datcompare.multi.long <- melt(datcompare.multi, id.vars = "peakid", variable.name = "Sample", value.name = "chic")

dat.multi.long <- dplyr::left_join(datcompare.multi.long, datref.multi)

# remove outliers
thres <- 0.995
dat.multi.long.filt <- dat.multi.long %>%
  group_by(Sample) %>%
  filter(chic < quantile(chic, probs = thres) & chic > quantile(chic, probs = 1 - thres)) %>%
  filter(chic_ref < quantile(chic_ref, probs = thres) & chic_ref > quantile(chic_ref, probs = 1 - thres))

ggplot(dat.multi.long.filt, aes(x = chic_ref, y = chic)) + geom_point(alpha = 0.25) + 
  scale_x_log10() + scale_y_log10() + facet_wrap(~Sample) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.cors <- dat.multi.long.filt %>%
  group_by(Sample) %>%
  summarise(corr.pears = cor(chic, chic_ref, method = "pearson"),
            corr.spear = cor(chic, chic_ref, method = "spearman"),
            corr.pears.log2 = cor(log2(chic), log2(chic_ref), method = "pearson"),
            corr.spear.log2 = cor(log2(chic), log2(chic_ref), method = "spearman"))
print(dat.cors)
