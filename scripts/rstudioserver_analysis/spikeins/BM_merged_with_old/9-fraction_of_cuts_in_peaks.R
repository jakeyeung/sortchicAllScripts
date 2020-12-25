# Jake Yeung
# Date of Creation: 2020-11-12
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/9-fraction_of_cuts_in_peaks.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(topicmodels)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

# Load outputs celltype specific peaks ------------------------------------

jmark <- "H3K4me1"
hubprefix <- "/home/jyeung/hub_oudenaarden"
indir <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/BM_merge_first_and_second.same_annot_file/merged_bams.first_and_second_rounds/counts_in_peaks_vs_nonpeaks_vs_blacklist.BM-AllMerged3")

# Annotate cells ----------------------------------------------------------

indir.annot <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.spikeins_mouse.BMround2_umaps_and_ratios.Round1Round2")
inf.annot <- file.path(indir.annot, paste0("spikeins_mouse.BMround1and2_umaps_and_ratios.colfix.celltyping.2020-11-01.WithRelLevels.mark_", jmark, ".cell_cluster_tables.txt"))
dat.annot <- fread(inf.annot)

# Load CSV by ctypes ------------------------------------------------------

jchromos <- paste(c(seq(19), "X", "Y"), sep = "")


dat.merge.annot.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  ctypes <- unique(sapply(list.files(indir, pattern = paste0(".*", jmark, ".*genome.csv")), function(x){
    # BM_round1_round2_merged_H3K4me1_pDCs.HiddenDomains.Eryths_2500_H3K4me1.cuts_in_genome.csv
    xsplit <- strsplit(strsplit(x, "\\.")[[1]][[1]], "_")[[1]][[6]]
    return(xsplit)
  }))
  names(ctypes) <- ctypes
  print(ctypes)

  dat.merge.annot.lst.bymark <- lapply(ctypes, function(jctype){
    print(jctype)
    
    # paste0("BM_round1_round2_merged_H3K9me3_HSPCs.HiddenDomains.Lymphoid_2500_H3K9me3.cuts_in_peaks.csv")
    inf1 <- list.files(indir, pattern = paste0("BM_round1_round2_merged_", jmark, "_", jctype, ".HiddenDomains.", jctype, ".*.cuts_in_peaks.csv"), full.names = TRUE)
    inf2 <- list.files(indir, pattern = paste0("BM_round1_round2_merged_", jmark, "_", jctype, ".HiddenDomains.", jctype, ".*.cuts_in_genome.csv"), full.names = TRUE)
    
    # jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
    dat1 <- fread(inf1) %>%
      group_by(samp) %>%
      summarise(cuts_in_peak = sum(cuts_in_peak))
    
    dat2 <- fread(inf2) %>%
      filter(chromosome %in% jchromos) %>%
      group_by(samp) %>%
      summarise(cuts_total = sum(cuts_in_chromo))
    
    spikeinchromo <- "J02459.1"
    spikeins <- fread(inf2) %>%
      filter(chromosome == spikeinchromo) %>%
      group_by(samp) %>%
      summarise(spikein_cuts = sum(cuts_in_blacklist))
    
    # NA spikeincuts if round1, nice! 
    dat.merge <- dat1 %>%
      left_join(., dat2, by = "samp")  %>%
      left_join(., spikeins, by = "samp") %>%
      ungroup() %>%
      mutate(ctype = jctype)
    
    dat.merge.annot <- left_join(dat.merge, dat.annot, by = c("samp" = "cell")) %>%
      rowwise() %>%
      mutate(mark = jmark,
             plate = ClipLast(samp, jsep = "_"))
    return(dat.merge.annot)
  }) %>%
    bind_rows()
}) 

mlst.total <- lapply(jmarks, function(jmark){
  print(jmark)
  m <- ggplot(dat.merge.annot.lst[[jmark]] %>% filter(!is.na(spikein_cuts)), aes(x = ctype, y = log2(cuts_total / spikein_cuts))) + 
    geom_boxplot() +
    geom_point() + 
    theme_bw() + 
    ggtitle(jmark) +
    facet_wrap(~plate) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
})

mlst <- lapply(jmarks, function(jmark){
  print(jmark)
  m <- ggplot(dat.merge.annot.lst[[jmark]], aes(x = ctype, y = log2(cuts_in_peak / (cuts_total - cuts_in_peak)))) + 
    geom_boxplot() +
    geom_point() + 
    theme_bw() + 
    ggtitle(jmark) +
    facet_wrap(~plate) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
})

mlst.spikeins <- lapply(jmarks, function(jmark){
  print(jmark)
  m <- ggplot(dat.merge.annot.lst[[jmark]] %>% filter(!is.na(spikein_cuts)), aes(x = ctype, y = log2(cuts_in_peak / spikein_cuts))) + 
    geom_boxplot() +
    geom_point() + 
    theme_bw() + 
    ggtitle(jmark) +
    facet_wrap(~plate) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  print(m)
})


# JFuncs::multiplot(mlst[[1]], mlst[[2]], mlst[[3]], mlst[[4]], cols = 2) 


# Get signal at each peak?  -----------------------------------------------




