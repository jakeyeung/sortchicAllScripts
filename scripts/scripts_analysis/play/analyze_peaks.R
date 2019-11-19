# Jake Yeung
# Date of Creation: 2019-03-13
# File: ~/projects/scchic/scripts/scripts_analysis/play/analyze_peaks.R
# Analyze Peaks for H3K27me3 and H3K9me3 make sure no label swap

library(cowplot)

# load data  --------------------------------------------------------------

# jmark <- "H3K27me3"
jmark <- "H3K9me3"

jmarks <- list("H3K9me3", "H3K27me3")
names(jmarks) <- jmarks

infs <- lapply(jmarks, function(jmark) paste0("/Users/yeung/data/scchic/beds/BM_", jmark, "_merged.1000.cutoff_analysis.blacklistfilt.annot.bed"))

dats <- lapply(jmarks, function(jmark){
  inf <- infs[[jmark]]
  dat <- data.table::fread(inf)
  colnames(dat) <- c("chromo", "start", "end", "gene", "genedist")
  dat <- dat %>%
    mutate(size = end - start)
  dat$mark <- jmark
  return(dat)
})

dats <- bind_rows(dats)

ggplot(dats, aes(x = size, fill = mark)) + geom_density(alpha = 0.5)

ggplot(dats, aes(x = abs(genedist), fill = mark)) + geom_density(alpha = 0.5) 
