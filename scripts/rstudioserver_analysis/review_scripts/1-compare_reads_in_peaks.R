# Jake Yeung
# Date of Creation: 2021-06-06
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/1-compare_reads_in_peaks.R
# 

rm(list=ls()) 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Constants  ------------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"
cnames <- c("fname1", "fname2", "nlines1", "nlines2")


# Filter by metadata good cells K562  -------------------------------------


inf.meta <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.K562_clean_by_G1filt.fracnonzerofilt/K562_cleaned_by_G1filt.H3K27me3.nonzerofracthres_2.bsize_50000.txt")
dat.meta <- fread(inf.meta)
good.cells <- dat.meta$cell


# Kaya-Okur ---------------------------------------------------------------


indir <- file.path(hubprefix, "jyeung/data/scChiC/public_data/K562_Kaya-Okur/H3K27me3/count_before_after_filtering_ENCODE_filt")
indir.ko.hdfilt <- file.path(hubprefix, "jyeung/data/scChiC/public_data/K562_Kaya-Okur/H3K27me3/count_before_after_filtering_hiddendomains_filt")

infs <- list.files(indir, pattern = "*.txt", full.names = TRUE)
infs.ko.hdfilt <- list.files(indir.ko.hdfilt, pattern = "*.txt", full.names = TRUE)

dat.ko.encodefilt <- lapply(infs, function(inf){
  fread(inf, header = FALSE)
}) %>%
  bind_rows()
colnames(dat.ko.encodefilt) <- cnames

dat.ko.encodefilt <- dat.ko.encodefilt %>%
  rowwise() %>%
  mutate(frac.reads.in.peaks = nlines2 / nlines1)

outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/K562_KayaOkur_vs_Zeller_H3K27me3.", Sys.Date(), ".pdf")
pdf(file = outpdf, useDingbats = FALSE)

ggplot(dat.ko.encodefilt, aes(x = frac.reads.in.peaks)) + 
  geom_density() 
dat.ko.encodefilt$samp <- sapply(dat.ko.encodefilt$fname1, function(x) basename(x))
 

# HD filt  ----------------------------------------------------------------


dat.ko.hdfilt <- lapply(infs.ko.hdfilt, function(inf){
  fread(inf, header = FALSE)
}) %>%
  bind_rows() 
colnames(dat.ko.hdfilt) <- cnames

dat.ko.hdfilt <- dat.ko.hdfilt %>%
  rowwise() %>%
  mutate(frac.reads.in.peaks = nlines2 / nlines1)
dat.ko.hdfilt$samp <- sapply(dat.ko.hdfilt$fname1, function(x) basename(x))

ggplot(dat.ko.hdfilt, aes(x = frac.reads.in.peaks)) + 
  geom_density() 


# Zeller  -----------------------------------------------------------------

inmain.z.encodefilt <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged")
indir.z.encodefilt <- file.path(inmain.z.encodefilt, "counts_in_peaks_vs_nonpeaks_vs_blacklist.faster.clean3.blfix.ChIP_seq_peaks")
# infs.z.encodefilt <- list.files(indir.z.encodefilt, pattern = "*reads_in_peaks.ChIP_seq_peaks.csv", full.names = TRUE)
inf.z.encodefilt.peaks <- file.path(indir.z.encodefilt, paste0("K562_AllMerged_H3K27me3.merged.sorted.tagged.bam.cuts_in_peaks.ChIP_seq_peaks.csv"))
inf.z.encodefilt.total <- file.path(indir.z.encodefilt, paste0("K562_AllMerged_H3K27me3.merged.sorted.tagged.bam.cuts_in_genome.ChIP_seq_peaks.csv"))

dat.z.encodefilt.peaks <- fread(inf.z.encodefilt.peaks) %>%
  group_by(samp) %>%
  summarise(cuts_in_peak = sum(cuts_in_peak)) %>%
  ungroup() %>%
  filter(samp %in% good.cells)

dat.z.encodefilt.total <- fread(inf.z.encodefilt.total) %>%
  group_by(samp) %>%
  summarise(cuts_in_chromo = sum(cuts_in_chromo)) %>%
  ungroup() %>%
  filter(samp %in% good.cells)

dat.z.encodefilt <- left_join(dat.z.encodefilt.peaks, dat.z.encodefilt.total, by = "samp") %>%
  rowwise() %>%
  mutate(frac.reads.in.peaks = cuts_in_peak / cuts_in_chromo)

ggplot(dat.z.encodefilt, aes(x = cuts_in_peak / cuts_in_chromo)) + 
  geom_density() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


# Compare Kaya vs Zeller  -------------------------------------------------

dat.compare <- rbind(subset(dat.ko.encodefilt, select = c(samp, frac.reads.in.peaks)) %>%
                           mutate(set = "KayaOkur"), 
                         subset(dat.z.encodefilt, select = c(samp, frac.reads.in.peaks)) %>%
                           mutate(set = "Zeller"))

ggplot(dat.compare, aes(x = frac.reads.in.peaks, fill = set)) + 
  geom_density(alpha = 0.5) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("Fraction of reads: ChIP-seq K562 ENCODE peaks")


# Compmare Kaya vs Zeller: HD filt  ---------------------------------------

inf.z.hdfilt.peaks <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/counts_in_peaks_vs_nonpeaks_vs_blacklist.faster.clean3.blfix/K562_AllMerged_H3K27me3.merged.sorted.tagged.bam.cuts_in_peaks.csv"))
inf.z.hdfilt.total <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_data_spikeins/K562_tech_rep_merged/counts_in_peaks_vs_nonpeaks_vs_blacklist.faster.clean3.blfix/K562_AllMerged_H3K27me3.merged.sorted.tagged.bam.cuts_in_genome.csv"))

dat.z.hdfilt.peaks <- fread(inf.z.hdfilt.peaks) %>%
  group_by(samp) %>%
  summarise(cuts_in_peak = sum(cuts_in_peak)) %>%
  ungroup() %>%
  filter(samp %in% good.cells)

dat.z.hdfilt.total <- fread(inf.z.hdfilt.total) %>%
  group_by(samp) %>%
  summarise(cuts_in_chromo = sum(cuts_in_chromo)) %>%
  ungroup() %>%
  filter(samp %in% good.cells)

dat.z.hdfilt <- left_join(dat.z.hdfilt.peaks, dat.z.hdfilt.total, by = "samp") %>%
  rowwise() %>%
  mutate(frac.reads.in.peaks = cuts_in_peak / cuts_in_chromo)


dat.compare.hdfilt <- rbind(subset(dat.ko.hdfilt, select = c(samp, frac.reads.in.peaks)) %>%
                           mutate(set = "KayaOkur"), 
                         subset(dat.z.hdfilt, select = c(samp, frac.reads.in.peaks)) %>%
                           mutate(set = "Zeller"))

ggplot(dat.compare.hdfilt, aes(x = frac.reads.in.peaks, fill = set)) + 
  geom_density(alpha = 0.5) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("Fraction of reads: HiddenDomains pseudobulk respective datasets")


# Compare sensitivities  --------------------------------------------------


ggplot(dat.z.hdfilt, aes(x = log2(cuts_in_chromo))) + 
  geom_density() +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(dat.ko.hdfilt, aes(x = log2(nlines1))) +
  geom_density() +
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


dat.hdfilt.merge <- rbind(dat.ko.hdfilt %>% 
                            dplyr::select(samp, nlines1, frac.reads.in.peaks) %>%
                            dplyr::rename(cuts_in_chromo = nlines1) %>%
                            dplyr::mutate(set = "KayaOkur"), 
                          dat.z.hdfilt %>%
                          dplyr::select(samp, cuts_in_chromo, frac.reads.in.peaks) %>%
                          dplyr::mutate(set = "Zeller"))

ggplot(dat.hdfilt.merge, aes(x = frac.reads.in.peaks, fill = set)) + 
  geom_density(alpha = 0.5) + 
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  xlab("Fraction of reads in peaks")

ggplot(dat.hdfilt.merge, aes(y = frac.reads.in.peaks, x = set, fill = set)) + 
  geom_boxplot(alpha = 0.5) + 
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ylab("Fraction of reads in peaks")

ggplot(dat.hdfilt.merge, aes(y = log2(frac.reads.in.peaks), x = set, fill = set)) + 
  geom_boxplot(alpha = 0.5) + 
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")

ggplot(dat.hdfilt.merge, aes(x = log2(cuts_in_chromo), fill = set)) + 
  geom_density(alpha = 0.5) + 
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  xlab("log2(unique reads)")

ggplot(dat.hdfilt.merge, aes(y = log2(cuts_in_chromo), x = frac.reads.in.peaks, color = set)) + 
  geom_point(alpha = 0.5) + 
  theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  ylab("log2(unique reads)") + 
  xlab("Fraction reads in peaks")

dev.off()
