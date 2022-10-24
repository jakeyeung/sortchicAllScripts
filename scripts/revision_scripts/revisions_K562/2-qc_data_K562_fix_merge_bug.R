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

outpdf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_cluster/K562_new/plots/qc_plots_K562_with_frip.", Sys.Date(), ".pdf")
outrds <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_cluster/K562_new/plots/qc_plots_K562_with_frip.", Sys.Date(), ".rds")


# Load QC data  -----------------------------------------------------------

jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks

inmain <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_cluster/K562_new/K562_match_marks"

dat.meta.lst <- lapply(jmarks, function(jmark){
  inf <- file.path(inmain, paste0("K562_", jmark), paste0("qc_metadata_new_only.K562_", jmark, ".0.8_0.5_3000.2022-03-11.txt"))
  print(inf)
  fread(inf) %>%
    mutate(mark = jmark)
})

m.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.meta.lst[[jmark]] %>%
    arrange(is.good) 
  ggplot(jdat, aes(x = log10(total.count.from.mat), y = TA.frac, color = is.good)) + 
    geom_vline(xintercept = log10(3000), linetype = "dashed", size = 0.25, alpha = 1, color = "grey55") + 
    geom_hline(yintercept = 0.5, linetype = "dashed", size = 0.25, alpha = 1, color = "grey55") + 
    geom_point(alpha = 0.2, size = 0.5) + 
    theme_bw() + 
    ggtitle(jmark) + 
    scale_y_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits = c(0, 1)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})

pdf(outpdf, useDingbats = FALSE)

print(JFuncs::multiplot(m.lst[[1]], m.lst[[2]], m.lst[[4]], m.lst[[3]], cols = 4))


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
print(m.plate.lst)

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
print(m.sum.lst)




# 
# Get fraction of cuts in peaks from hiddendomains  -----------------------

# indir.frips <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_cluster/K562_new/counts_in_peaks_vs_nonpeaks_vs_blacklist.faster.clean3.blfix"
indir.frips <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_cluster/K562_new/counts_in_peaks_vs_nonpeaks_vs_blacklist.faster.clean3.blfix.with_bad_cells.nochr"


dat.merge.lst <- lapply(jmarks, function(jmark){
  inf.genome <- file.path(indir.frips, paste0("K562_", jmark, "_plates_merged.tagged.bam.cuts_in_genome.csv"))
  print(inf.genome)
  dat.genome <- fread(inf.genome) %>%
    group_by(samp) %>%
    summarise(total_cuts = sum(cuts_in_chromo))
  inf.peaks <- file.path(indir.frips, paste0("K562_", jmark, "_plates_merged.tagged.bam.cuts_in_peaks.csv"))
  print(inf.peaks)
  dat.peaks <- fread(inf.peaks) %>%
    group_by(samp) %>%
    summarise(total_cuts_in_peaks = sum(cuts_in_peak))
  dat.merged <- left_join(dat.genome, dat.peaks) %>%
    left_join(., dat.meta.lst[[jmark]])
})

m.lst.frips <- lapply(jmarks, function(jmark){
  print(head(dat.merge.lst[[jmark]]))
  ggplot(dat.merge.lst[[jmark]], aes(x = total_cuts_in_peaks / total_cuts, fill = is.good)) + 
    geom_density() + 
    theme_bw() + 
    ggtitle(jmark) + 
    scale_x_continuous(breaks = c(0, 0.5, 1), labels = c(0, 0.5, 1), limits = c(0, 1)) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
})

library(JFuncs)
multiplot(m.lst.frips[[1]], m.lst.frips[[2]], m.lst.frips[[4]], m.lst.frips[[3]], cols = 4)

saveRDS(dat.merge.lst, file = outrds)

dev.off()



# Calculate fraction of cells that are "good" ------------------------------

dat.meta.merge <- do.call(rbind, dat.meta.lst)

dat.is.good.frac.mark <- dat.meta.merge %>%
  group_by(is.good, mark) %>%
  summarise(ncells = length(samp)) %>%
  group_by(mark) %>%
  mutate(nfrac = ncells / sum(ncells))

dat.is.good.frac <- dat.meta.merge %>%
  group_by(is.good) %>%
  summarise(ncells = length(samp)) %>%
  ungroup() %>%
  mutate(nfrac = ncells / sum(ncells))

dat.meta.merge.good <- dat.meta.merge %>%
  filter(is.good)

plot(density(log10(dat.meta.merge.good$total.count.from.mat)))
mean(dat.meta.merge.good$total.count.from.mat)
mean(dat.meta.merge.good$total.count)
median(dat.meta.merge.good$total.count.from.mat)
median(dat.meta.merge.good$total.count)

median(dat.meta.merge.good$frac.count.in.peak)
