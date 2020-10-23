# Jake Yeung
# Date of Creation: 2020-08-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged_for_MARA/compare_MARA_ZF_and_BM.R
# 


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(ggrepel)

# Load maras --------------------------------------------------------------

# jmark <- "H3K4me1"
# jmark <- "H3K27me3"
jmark <- "H3K4me3"
  
mdir.BM <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/", jmark, "/mara_output/BM_", jmark, ".BM_AllMerged.VarCorrection.Poisson.GLMPCA-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, "/BM_", jmark, ".BM_AllMerged.VarCorrection.Poisson.GLMPCA")
assertthat::assert_that(dir.exists(mdir.BM))
mara.out.BM <- LoadMARA(mdir = mdir.BM, make.cnames = FALSE)

mdir.ZF <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_ZF-AllMerged2_Peaks_1000/", jmark, "/mara_output/ZF_", jmark, ".ZF_AllMerged.VarCorrection.Poisson.GLMPCA.txt-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, "/ZF_", jmark, ".ZF_AllMerged.VarCorrection.Poisson.GLMPCA.txt")
assertthat::assert_that(dir.exists(mdir.ZF))
mara.out.ZF <- LoadMARA(mdir = mdir.ZF, make.cnames = FALSE)


dat.zscores.BM <- data.frame(motif = mara.out.BM$zscore$motif, zscore.BM = mara.out.BM$zscore$zscore, stringsAsFactors = FALSE)
dat.zscores.ZF <- data.frame(motif = mara.out.ZF$zscore$motif, zscore.ZF = mara.out.ZF$zscore$zscore, stringsAsFactors = FALSE)

dat.zscores <- left_join(dat.zscores.BM, dat.zscores.ZF, by = "motif") %>%
  rowwise() %>%
  mutate(motif.lab = ifelse(zscore.BM > 0.5 | zscore.ZF > 0.5, as.character(motif), NA))

ggplot(dat.zscores, aes(x = zscore.BM, y = zscore.ZF, label = motif.lab)) + 
  geom_point()  + 
  geom_text_repel() 




# Copare two marks --------------------------------------------------------

jmark1 <- "H3K4me1"
jmark2 <- "H3K4me3"

mdir.ZF1 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_ZF-AllMerged2_Peaks_1000/", jmark1, "/mara_output/ZF_", jmark1, ".ZF_AllMerged.VarCorrection.Poisson.GLMPCA.txt-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark1, "/ZF_", jmark1, ".ZF_AllMerged.VarCorrection.Poisson.GLMPCA.txt")
assertthat::assert_that(dir.exists(mdir.ZF1))
mara.out.ZF1 <- LoadMARA(mdir = mdir.ZF1, make.cnames = FALSE)

mdir.ZF2 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_ZF-AllMerged2_Peaks_1000/", jmark2, "/mara_output/ZF_", jmark2, ".ZF_AllMerged.VarCorrection.Poisson.GLMPCA.txt-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark2, "/ZF_", jmark2, ".ZF_AllMerged.VarCorrection.Poisson.GLMPCA.txt")
assertthat::assert_that(dir.exists(mdir.ZF2))
mara.out.ZF2 <- LoadMARA(mdir = mdir.ZF2, make.cnames = FALSE)

dat.zscores.ZF1 <- data.frame(motif = mara.out.ZF1$zscore$motif, zscore.ZF1 = mara.out.ZF1$zscore$zscore, stringsAsFactors = FALSE)
dat.zscores.ZF2 <- data.frame(motif = mara.out.ZF2$zscore$motif, zscore.ZF2 = mara.out.ZF2$zscore$zscore, stringsAsFactors = FALSE)

dat.zscores.ZF.marks <- left_join(dat.zscores.ZF1, dat.zscores.ZF2, by = "motif") %>%
  rowwise() %>%
  mutate(motif.lab = ifelse(zscore.ZF1 > 0.5 | zscore.ZF2 > 0.5, as.character(motif), NA))

ggplot(dat.zscores.ZF.marks, aes(x = zscore.ZF1, y = zscore.ZF2, label = motif.lab)) + 
  geom_point()  + 
  geom_text_repel() 



# Compare BM --------------------------------------------------------------


jmark1 <- "H3K4me1"
jmark2 <- "H3K4me3"

# mdir.BM1 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_ZF-AllMerged2_Peaks_1000/", jmark1, "/mara_output/ZF_", jmark1, ".ZF_AllMerged.VarCorrection.Poisson.GLMPCA.txt-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark1, "/ZF_", jmark1, ".ZF_AllMerged.VarCorrection.Poisson.GLMPCA.txt")
mdir.BM1 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/", jmark1, "/mara_output/BM_", jmark1, ".BM_AllMerged.VarCorrection.Poisson.GLMPCA-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark1, "/BM_", jmark1, ".BM_AllMerged.VarCorrection.Poisson.GLMPCA")
assertthat::assert_that(dir.exists(mdir.BM1))
mara.out.BM1 <- LoadMARA(mdir = mdir.BM1, make.cnames = FALSE)

mdir.BM2 <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/mara_analysis_BM-AllMerged_Peaks_1000/", jmark2, "/mara_output/BM_", jmark2, ".BM_AllMerged.VarCorrection.Poisson.GLMPCA-hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark2, "/BM_", jmark2, ".BM_AllMerged.VarCorrection.Poisson.GLMPCA")
assertthat::assert_that(dir.exists(mdir.BM2))
mara.out.BM2 <- LoadMARA(mdir = mdir.BM2, make.cnames = FALSE)


dat.zscores.BM1 <- data.frame(motif = mara.out.BM1$zscore$motif, zscore.BM1 = mara.out.BM1$zscore$zscore, stringsAsFactors = FALSE)
dat.zscores.BM2 <- data.frame(motif = mara.out.BM2$zscore$motif, zscore.BM2 = mara.out.BM2$zscore$zscore, stringsAsFactors = FALSE)

dat.zscores.BM.marks <- left_join(dat.zscores.BM1, dat.zscores.BM2, by = "motif") %>%
  rowwise() %>%
  mutate(motif.lab = ifelse(zscore.BM1 > 0.5 | zscore.BM2 > 0.5, as.character(motif), NA))

ggplot(dat.zscores.BM.marks, aes(x = zscore.BM1, y = zscore.BM2, label = motif.lab)) + 
  geom_point()  + 
  geom_text_repel() 



