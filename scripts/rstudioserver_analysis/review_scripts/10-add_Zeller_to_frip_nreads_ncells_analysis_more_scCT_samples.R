# Jake Yeung
# Date of Creation: 2021-06-22
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/10-add_Zeller_to_frip_nreads_ncells_analysis_more_scCT_samples.R
# 


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)



# Load zeller  ------------------------------------------------------------

# inf.barcodes <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/unique_reads_many_studies.2021-06-13.Zeller_dedup_fixed.merged_with_cellspec_norm.rds"
# dat.barcodes <- readRDS(inf.barcodes)

inf.zeller <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/dat_hdfilt_merge.KayaOkur_vs_Zeller.H3K27me3.2021-06-07.rds"
dat.zeller <- readRDS(inf.zeller) %>%
  # filter(set == "Zeller") %>%
  dplyr::rename("cell" = "samp", 
                allcounts = "cuts_in_chromo", 
                jset = "set") %>%
  rowwise() %>%
  mutate(paper = jset, 
         peakcounts = allcounts * frac.reads.in.peaks) %>%
  dplyr::select(-frac.reads.in.peaks)
  


# Load bart ---------------------------------------------------------------

jdate <- "2021-06-22"
indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data"
inname <- paste0("frip_peakcounts_allcounts.all_scCT.", jdate, ".rds")
inrds <- file.path(indir, inname)
dat.bart <- readRDS(inrds)

dat.bartzeller <- bind_rows(dat.zeller, dat.bart)

dat.bartzeller.filt <- dat.bartzeller %>%
  filter( (paper == "grosselin" & allcounts > 1000) | (paper != "grosselin") ) %>%
  filter(paper != "kaya_okur")

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/reanalysis_Bartosovic"
pdfname <- paste0("bartosovic_reanalysis_with_Zeller_more_scCT.", Sys.Date(), ".pdf")
outpdf <- file.path(outdir, pdfname)

pdf(outpdf, useDingbats = FALSE)

ggplot(dat.bartzeller.filt, aes(x = jset, y = peakcounts / allcounts, fill = paper)) + 
  geom_boxplot() + 
  theme_bw() + 
  ylab("Fraction of Unique Reads in Peaks") + 
  xlab("") + 
  scale_fill_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.bartzeller.filt, aes(x = jset, y = peakcounts / allcounts, fill = paper)) + 
  geom_boxplot() + 
  theme_bw() + 
  scale_y_log10() + 
  ylab("Fraction of Unique Reads in Peaks") + 
  xlab("") + 
  scale_fill_manual(values = cbPalette) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# get nreads 
ggplot(dat.bartzeller.filt, aes(x = jset, y = allcounts, fill = paper)) + 
  scale_y_log10() + 
  geom_boxplot() + 
  theme_bw() + 
  ylab("Unique Reads") + 
  xlab("") + 
  scale_fill_manual(values = cbPalette) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# get ncells
dat.bartzeller.sum.filt <- dat.bartzeller.filt %>%
  group_by(jset, paper) %>%
  summarise(ncells = length(cell))

ggplot(dat.bartzeller.sum.filt, aes(x = jset, y = log2(ncells), fill = paper)) + 
  geom_col() + 
  theme_bw() + 
  scale_fill_manual(values = cbPalette) +
  xlab("") + 
  scale_y_continuous(breaks = seq(12)) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.bartzeller.sum.filt, aes(x = jset, y = ncells, fill = paper)) + 
  geom_col() + 
  xlab("") + 
  theme_bw() + 
  scale_fill_manual(values = cbPalette) +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

dev.off()

outrdsmerged <- file.path(outdir, paste0("Bartosovic_reanalysis_many_datasets_frip_with_Zeller_more_scCT.", Sys.Date(), ".rds"))
saveRDS(dat.bartzeller.filt, file = outrdsmerged)

