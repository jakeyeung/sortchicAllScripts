# Jake Yeung
# Date of Creation: 2020-08-05
# File: ~/projects/scchic/scripts/rstudioserver_analysis/dinuc_analysis/init_dinuc_analysis3_tlen147.R
# 


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

DoFFT <- function(y){
  yt <- fft(y)[1:(length(y) / 2)]
  xt <- (seq(yt) - 1) / length(y)
  fout <- data.frame(k = seq(yt) - 1, f = xt, y = yt, amp = Mod(yt), period = 1 / xt, stringsAsFactors = FALSE)
  return(fout)
}


# Granu -------------------------------------------------------------------

indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/dinuc_freq/tlen_220/freq_cuts_downstream"
# indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/dinuc_freq/tlen_190/freq_cuts_downstream"
inf <- file.path(indir, "H3K27me3-BM_AllMerged.HSPCs.sorted.tlenfilt.csv.dinucfreq.csv.gz")
# inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/dinuc_freq/tlen_147_reverse/freq_cuts_downstream/H3K27me3-BM_AllMerged.HSPCs.sorted.tlenfilt.csv.dinucfreq.csv.gz"

dat <- fread(inf, sep = "\t", header = FALSE)

dat.sum <- colMeans(dat)

# plot(dat.sum[50:150])
# plot(dat.sum[10:100])

plot(dat.sum, type = "o")
abline(v = 147)
abline(v = 74)
abline(v = 190)

plot(dat.sum, type = "o")
plot(dat.sum[20:175], type = "o")

# yspace <- dat.sum[1:150]
yspace <- dat.sum[150:200]

fft.out <- DoFFT(y = yspace)

ggplot(fft.out %>% filter(k > 0), aes(x = f, y = amp)) + geom_point() 

# 
# # Check frag sizes --------------------------------------------------------
# 
# 
# indir.frags <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data/ZellerRawData_B6_All_MergedByMarks_final.bams_by_cluster.MAPQ_40.2020-06-17/bamPEFragmentSize"
# infs.frags <- list.files(indir.frags, pattern = "*.txt", full.names = TRUE)
# 
# dat.lst <- lapply(infs.frags, function(inf){
#   dat.tmp <- fread(inf)
#   return(dat.tmp)
# }) %>%
#   bind_rows() %>%
#   group_by(Sample) %>% 
#   filter(Size > 0 & Size < 1000) %>%
#   mutate(Freq = Occurrences / sum(Occurrences))
# 
# ggplot(dat.lst, aes(x = Size, y = Freq)) + 
#   geom_col(width = 2, alpha = 0.5) + 
#   theme_bw() + 
#   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_wrap(~Sample) + 
#   geom_vline(xintercept = c(147, 190, 220), color = 'red')
# 
# 
# dat.bymark <- dat.lst %>%
#   rowwise() %>%
#   mutate(mark = strsplit(Sample, split = "-")[[1]][[1]]) %>%
#   group_by(mark, Size) %>%
#   summarise(Occurrences = sum(Occurrences)) %>%
#   group_by(mark) %>%
#   mutate(Freq = Occurrences / sum(Occurrences))
# 
# ggplot(dat.bymark, aes(x = Size, y = Freq)) + 
#   geom_point(alpha = 0.5) +
#   theme_bw() + 
#   theme(aspect.ratio=0.25, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_wrap(~mark, ncol = 1) + 
#   geom_vline(xintercept = c(147, 190, 220), color = 'blue')
# 
# 
# 
# 
# 
# 
# 
