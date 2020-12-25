# Jake Yeung
# Date of Creation: 2020-11-29
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/15-compare_H3K27me3_spikeins_reseq.R
# Compare H3K27me3 amounts by celltype before and after reseq


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)


hubprefix <- "/home/jyeung/hub_oudenaarden"
jchromos <- paste(c(seq(19), "X", "Y"), sep = "")
jmarks <- c("H3K27me3"); names(jmarks) <- jmarks

spikeinchromo <- "J02459.1"


# Load meta ---------------------------------------------------------------

inf.meta <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2.redo_after_varfilt_with_spikeins/cell_cluster_table_with_spikeins.H3K27me3.2020-11-18.dupfilt.txt"
dat.meta <- fread(inf.meta)

dat.annot <- subset(dat.meta, select = c(cell, cluster, Cluster, batch, plate))

# Load reseq spikeins --------------------------------------------------------------

# indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_H3K27me3_tech_rep_merged/raw_demultiplexed/tagged_bams/counts_in_peaks_vs_nonpeaks_vs_blacklist.BM-AllMerged3.peaks_merged.reseq")
indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/split_by_cluster.MAPQ_40.H3K27me3reseq/H3K27me3/merged_by_ctype/counts_in_peaks_vs_nonpeaks_vs_blacklist")


dat.merge.annot.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  ctypes <- unique(sapply(list.files(indir, pattern = paste0(".*", jmark, ".*genome.csv")), function(x){
    # PZ-BM-rep3-H3K27me3-platemerged.Eryths.HiddenDomains.Eryths_2500_H3K27me3.cuts_in_peaks.csv
    # xsplit <- strsplit(strsplit(x, "\\.")[[1]][[1]], "-")[[1]][[5]]
    xsplit <- strsplit(x, "\\.")[[1]][[2]]
    return(xsplit)
  }))
  names(ctypes) <- ctypes
  print(ctypes)
  
  dat.merge.annot.lst.bymark <- lapply(ctypes, function(jctype){
    print(jctype)
    
    # paste0("BM_round1_round2_merged_H3K9me3_HSPCs.HiddenDomains.Lymphoid_2500_H3K9me3.cuts_in_peaks.csv")
    # # PZ-BM-rep3-H3K27me3-8.sorted.tagged.HiddenDomains.__H3K27me3.cuts_in_genome.csv 
    # PZ-BM-rep3-H3K27me3-platemerged.Granulocytes.HiddenDomains.Granulocytes_2500_H3K27me3.cuts_in_peaks.csv 
    inf1 <- list.files(indir, pattern = paste0("PZ-BM-rep3-", jmark, "-platemerged.", jctype, ".HiddenDomains.", jctype, ".*.cuts_in_peaks.csv"), full.names = TRUE)  # * because 2500 or 5000 is possible
    inf2 <- list.files(indir, pattern = paste0("PZ-BM-rep3-", jmark, "-platemerged.", jctype, ".HiddenDomains.", jctype, ".*.cuts_in_genome.csv"), full.names = TRUE)  # * because 2500 or 5000 is possible
    
    # jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
    dat1 <- fread(inf1) %>%
      group_by(samp) %>%
      summarise(cuts_in_peak2 = sum(cuts_in_peak))
    
    dat2 <- fread(inf2) %>%
      filter(chromosome %in% jchromos) %>%
      group_by(samp) %>%
      summarise(cuts_total2 = sum(cuts_in_chromo))
    
    spikeins <- fread(inf2) %>%
      filter(chromosome == spikeinchromo) %>%
      group_by(samp) %>%
      summarise(spikein_cuts2 = sum(cuts_in_blacklist))
    
    # NA spikeincuts if round1, nice! 
    dat.merge <- dat1 %>%
      left_join(., dat2, by = "samp")  %>%
      left_join(., spikeins, by = "samp") %>%
      ungroup() %>%
      mutate(ctype = jctype)
    
    dat.merge.annot <- left_join(dat.merge, dat.annot, by = c("samp" = "cell"))
    return(dat.merge.annot)
  }) %>%
    bind_rows()
}) 



head(dat.meta)


# Compare numbers before and after  ---------------------------------------


dat.meta.merge <- left_join(dat.merge.annot.lst[[jmarks[[1]]]], subset(dat.meta, select = c(cell, cuts_in_peak, cuts_total, spikein_cuts)), by = c("samp" = "cell"))

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

ggplot(dat.meta.merge, aes(x = cuts_in_peak, y = cuts_in_peak2, color = plate)) + 
  geom_point()  + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  facet_wrap(~plate) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_abline(slope = 1, intercept = 0)

ggplot(dat.meta.merge, aes(x = cuts_in_peak, y = cuts_in_peak2, color = plate)) + 
  geom_point()  + 
  scale_x_log10() + 
  scale_y_log10() + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  facet_wrap(~plate) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_abline(slope = 1, intercept = 0)

ggplot(dat.meta.merge, aes(x = cuts_total, y = cuts_total2, color = plate)) + 
  geom_point()  + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  facet_wrap(~plate) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_abline(slope = 1, intercept = 0)

ggplot(dat.meta.merge, aes(x = cuts_total, y = cuts_total2, color = plate)) + 
  geom_point()  + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  facet_wrap(~plate) + 
  scale_x_log10() + 
  scale_y_log10() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_abline(slope = 1, intercept = 0)

ggplot(dat.meta.merge, aes(x = spikein_cuts, y = spikein_cuts2, color = plate)) + 
  geom_point()  + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  facet_wrap(~plate) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_abline(slope = 1, intercept = 0)

ggplot(dat.meta.merge, aes(x = spikein_cuts, y = spikein_cuts2, color = plate)) + 
  geom_point()  + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  facet_wrap(~plate) + 
  scale_x_log10() + 
  scale_y_log10() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_abline(slope = 1, intercept = 0)

ggplot(dat.meta.merge, aes(x = cuts_in_peak / spikein_cuts, y = cuts_in_peak2 / spikein_cuts2, color = plate)) + 
  geom_point()  + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  facet_wrap(~plate) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_abline(slope = 1, intercept = 0)


ggplot(dat.meta.merge, aes(x = log2(cuts_in_peak / spikein_cuts), y = log2(cuts_in_peak2 / spikein_cuts2), color = plate)) + 
  geom_point()  + 
  scale_color_manual(values = cbPalette) + 
  theme_bw() + 
  facet_wrap(~plate) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  geom_abline(slope = 1, intercept = 0)


