# Jake Yeung
# Date of Creation: 2022-03-22
# File: ~/projects/scchic/scripts/revision_scripts/revisions_K562/2-qc_data_K562.R
# 

rm(list=ls())

  
library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# Load QC data  -----------------------------------------------------------

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks

inmain <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_cluster/K562_new/K562_match_marks"

dat.meta.lst <- lapply(jmarks, function(jmark){
  inf <- file.path(inmain, paste0("K562_", jmark), paste0("qc_metadata_new_only.K562_", jmark, ".0.8_0.5_3000.2022-03-11.txt"))
  fread(inf) %>%
    mutate(mark = jmark)
})

m.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.meta.lst[[jmark]] %>%
    arrange(is.good) 
  ggplot(jdat, aes(x = log10(total.count.from.mat), y = TA.frac, color = is.good)) + 
    geom_point(alpha = 0.2, size = 1) + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})


m.plate.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.meta.lst[[jmark]] %>%
    arrange(is.good) 
  ggplot(jdat, aes(x = log10(total.count.from.mat), y = TA.frac, color = is.good)) + 
    geom_point() + 
    facet_wrap(~plate) + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
m.plate.lst

m.sum.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.meta.lst[[jmark]] %>%
    arrange(is.good) 
  ggplot(jdat, aes(x = log10(total.count.from.mat))) + 
    geom_density() + 
    # facet_wrap(~plate) + 
    theme_bw() + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})
m.sum.lst


# Get fraction of cuts in peaks from hiddendomains  -----------------------

# indir.frips <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_cluster/K562_new/counts_in_peaks_vs_nonpeaks_vs_blacklist.faster.clean3.blfix"
indir.frips <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_cluster/K562_new/counts_in_peaks_vs_nonpeaks_vs_blacklist.faster.clean3.blfix"

# test
jmark <- "k9me3"

dat.genome.lst <- lapply(jmarks, function(jmark){
  inf.genome <- file.path(indir.frips, paste0("K562_", jmark, "_plates_merged.tagged.GoodCells.sorted.bam.cuts_in_genome.csv"))
  print(inf.genome)
  dat.genome <- fread(inf.genome) %>%
    group_by(chromosome, samp) %>%
    summarise(total_cuts = sum(cuts_in_chromo))
})


dat.peaks.lst <- lapply(jmarks, function(jmark){
  inf.peaks <- file.path(indir.frips, paste0("K562_", jmark, "_plates_merged.tagged.GoodCells.sorted.bam.cuts_in_peaks.csv"))
  print(inf.peaks)
  dat.peaks <- fread(inf.peaks) %>%
    group_by(chromosome, samp) %>%
    summarise(total_cuts_in_peaks = sum(cuts_in_peak))
})

plot(density(log10(dat.test.sum$total_cuts)), xlab = "log10(total counts)", main = "H3K9me3")


# Filter out meta good cells only  ----------------------------------------


metadirout <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_cluster/K562_new/K562_match_marks_meta_good_cells_only"
dir.create(metadirout)
for (jmark in jmarks){
  print(jmark)
  jtmp <- subset(dat.meta.lst[[jmark]], is.good)
  ftmp <- file.path(metadirout, paste0("metadata_", jmark, "_good_cells_only.", Sys.Date(), ".txt"))
  fwrite(jtmp, ftmp, sep = "\t")
}
