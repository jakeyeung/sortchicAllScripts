# Jake Yeung
# Date of Creation: 2021-06-16
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/7-check_Grosselin.R
# 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

hubprefix <- "/home/jyeung/hub_oudenaarden"

inf.tcell <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Grosselin_et_al_2019/GSE117309_RAW/GSM3290888_CountTable_Jurkat_scChIP_K27me3.txt.gz")
dat.tcell <- fread(inf.tcell)[, -1]

inf.bcell <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Grosselin_et_al_2019/GSE117309_RAW/GSM3290888_CountTable_Ramos_scChIP_K27me3.txt.gz")
dat.bcell <- fread(inf.bcell)[, -1]

dat.tcell.sum <- data.frame(cell = colnames(dat.tcell), nreads = colSums(dat.tcell), jset = "Grosselin_JurkatTcell", stringsAsFactors = FALSE)
dat.bcell.sum <- data.frame(cell = colnames(dat.bcell), nreads = colSums(dat.bcell), jset = "Grosselin_RamosBcell", stringsAsFactors = FALSE)
dat.grosselin.sum <- rbind(dat.tcell.sum, dat.bcell.sum)

ggplot(dat.grosselin.sum, aes(x = jset, y = nreads)) + 
  geom_point() + 
  geom_boxplot() + 
  scale_y_log10() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())



# Integrate ---------------------------------------------------------------


inf.ref <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/unique_reads_many_studies.2021-06-13.Zeller_dedup_fixed.rds"
dat.ref <- readRDS(inf.ref)

dat.merged <- rbind(dat.ref, dat.grosselin.sum)


jset2tech <- hash::hash(c("Zeller", "Ku_WBCs", "KayaOkur", "Grosselin_RamosBcell", "Grosselin_JurkatTcell", "Bart2021_Oli-neu", "Bart2021_mESC", "Bart2021_3T3", "Wu2021_hESC"), 
                        c("MNase", "MNase", "Tn5", "ChIP", "ChIP", "Tn5", "Tn5", "Tn5", "Tn5"))

dat.merged$tech <- sapply(as.character(dat.merged$jset), function(x) jset2tech[[x]][[1]])

cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# ggplot(dat.merged, aes(x = jset, y = nreads)) + 
ggplot(dat.merged, aes(x = forcats::fct_reorder(jset, nreads, median, .desc = TRUE), y = nreads, fill = tech)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cbPalette) + 
  scale_y_log10() + 
  xlab("") + ylab("Unique Reads") + 
  theme_bw(24) + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


inf.reftotal <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/post_submission/K562_public_data/unique_reads_many_studies.2021-06-13.Zeller_dedup_fixed.merged_with_cellspec_norm.rds"
dat.reftotal <- readRDS(inf.reftotal)

# normalize by average ntotal 
dat.normfacs <- dat.reftotal %>%
  group_by(jset) %>%
  summarise(medtotal = mean(ntotal),
            ncells = length(cell),
            normfac = medtotal / ncells)

dat.reftotal.normfacs <- left_join(dat.reftotal, dat.normfacs)

ggplot(dat.reftotal.normfacs, aes(x = forcats::fct_reorder(jset, nreads, median, .desc = TRUE), y = nreads)) + 
  geom_point(alpha = 0.5) +  
  geom_boxplot() + 
  scale_y_log10() + 
  theme_bw() + 
  xlab("") + ylab("Unique reads per cell / AverageLibrarySizePerCell") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

ggplot(dat.reftotal.normfacs, aes(x = forcats::fct_reorder(jset, nreads, median, .desc = TRUE), y = nreads / normfac)) + 
  geom_point(alpha = 0.5) +  
  geom_boxplot() + 
  scale_y_log10() + 
  theme_bw() + 
  xlab("") + ylab("Unique reads per cell / AverageLibrarySizePerCell") + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))





