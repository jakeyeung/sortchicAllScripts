# Jake Yeung
# Date of Creation: 2021-06-13
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/4-doublecheck_Wu.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

hubprefix <- "/home/jyeung/hub_oudenaarden"

inf.hesc <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Wu_et_al_2021/from_github/K27me3_stdcells.fragments.tsv.gz")
inf.h1de <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Wu_et_al_2021/from_github/K27me3_h1de.fragments.tsv.gz")

cnames <- c("chromo", "start", "end", "cell", "ntotal")
dat.hesc <- fread(inf.hesc, header = FALSE, col.names = cnames)
dat.h1de <- fread(inf.h1de, header = FALSE, col.names = cnames)

dat.hesc.sum <- dat.hesc %>%
  rowwise() %>%
  mutate(uniqCounts = 1) %>%
  group_by(cell) %>%
  summarise(ntotal = sum(ntotal), 
            nreads = sum(uniqCounts)) %>%
  mutate(jset = "Wu2021_hESCs")

dat.h1de.sum <- dat.h1de %>%
  rowwise() %>%
  mutate(uniqCounts = 1) %>%
  group_by(cell) %>%
  summarise(ntotal = sum(ntotal), 
            nreads = sum(uniqCounts)) %>%
  mutate(jset = "Wu2021_H1DE")

dat.wu <- rbind(dat.hesc.sum, dat.h1de.sum)

ggplot(dat.wu, aes(x = jset, y = log10(nreads))) + 
  geom_boxplot() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.wu, aes(x = ntotal, y = nreads, color = jset)) + 
  scale_x_log10() + 
  scale_y_log10() + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


inrds.merged <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/unique_reads_many_studies.2021-06-13.Zeller_dedup_fixed.merged_with_cellspec_norm.rds")
dat.all.merge <- readRDS(inrds.merged)

dat.all.noWu <- subset(dat.all.merge, jset != "Wu2021_hESC")

dat.all.merge.fix <- bind_rows(dat.all.noWu, dat.wu)

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.all.merge.fix, aes(x = ntotal, y = nreads, color = jset)) + 
  geom_point(alpha = 0.5) +  
  scale_x_log10() + 
  scale_y_log10() + 
  theme_bw() + 
  scale_color_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggplot(dat.all.merge.fix, aes(x = forcats::fct_reorder(jset, nreads, median, .desc = TRUE), y = nreads)) + 
  geom_point(alpha = 0.5) +  
  geom_boxplot() + 
  scale_y_log10() + 
  theme_bw() + 
  xlab("") + ylab("Unique reads") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# normalize by average ntotal 
dat.normfacs <- dat.all.merge.fix %>%
  group_by(jset) %>%
  summarise(normfac = median(ntotal))

dat.all.merge.fix.norm <- left_join(dat.all.merge.fix, dat.normfacs)
  
ggplot(dat.all.merge.fix.norm, aes(x = forcats::fct_reorder(jset, nreads, median, .desc = TRUE), y = nreads / normfac)) + 
  geom_point(alpha = 0.5) +  
  geom_boxplot() + 
  scale_y_log10() + 
  theme_bw() + 
  xlab("") + ylab("Unique reads per cell / AverageLibrarySizePerCell") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


ggplot(dat.all.merge.fix.norm, aes(x = forcats::fct_reorder(jset, nreads, median, .desc = TRUE), y = nreads / ntotal)) + 
  geom_point(alpha = 0.5) +  
  geom_boxplot() + 
  scale_y_log10() + 
  theme_bw() + 
  xlab("") + ylab("Unique reads, normalize by cell-specific library size") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))



# Integrate with other data  ----------------------------------------------


