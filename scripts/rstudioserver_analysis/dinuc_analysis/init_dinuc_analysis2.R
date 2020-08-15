# Jake Yeung
# Date of Creation: 2020-08-05
# File: ~/projects/scchic/scripts/rstudioserver_analysis/dinuc_analysis/init_dinuc_analysis.R
# 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


# Granu -------------------------------------------------------------------

# inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/dinuc_freq/again2/freq_cuts_downstream2/H3K4me3-BM_AllMerged.Bcells.sorted.cleaned.csv.dinucfreq.csv"

inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/dinuc_freq/again2/freq_cuts_downstream2/H3K27me3-BM_AllMerged.HSPCs.sorted.cleaned.csv.dinucfreq.csv"
dat <- fread(inf, sep = "\t", header = FALSE)

dat.sum <- colMeans(dat)

# plot(dat.sum[50:150])
plot(dat.sum)


# Check frag sizes --------------------------------------------------------


indir.frags <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/bamPEFragmentSize"
infs.frags <- list.files(indir.frags, pattern = "*.txt", full.names = TRUE)

dat.lst <- lapply(infs.frags, function(inf){
  dat.tmp <- fread(inf)
  return(dat.tmp)
}) %>%
  bind_rows() %>%
  group_by(Sample) %>% 
  filter(Size > 0 & Size < 1000) %>%
  mutate(Freq = Occurrences / sum(Occurrences))

ggplot(dat.lst, aes(x = Size, y = Freq)) + 
  geom_col(width = 2, alpha = 0.5) + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~Sample) + 
  geom_vline(xintercept = c(147, 190, 220), color = 'red')


dat.bymark <- dat.lst %>%
  rowwise() %>%
  mutate(mark = strsplit(Sample, split = "-")[[1]][[1]]) %>%
  group_by(mark, Size) %>%
  summarise(Occurrences = sum(Occurrences)) %>%
  group_by(mark) %>%
  mutate(Freq = Occurrences / sum(Occurrences))

ggplot(dat.bymark, aes(x = Size, y = Freq)) + 
  geom_point(alpha = 0.5) +
  theme_bw() + 
  theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~mark, ncol = 1) + 
  geom_vline(xintercept = c(147, 190, 220), color = 'blue')







