# Jake Yeung
# Date of Creation: 2020-08-27
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged_for_MARA/remove_large_peaks_write_clean_N_and_E.BM.R
# 


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)

library(mixtools)

# Explore peaks  ----------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

# jmark <- "H3K4me1"
# jmark <- "H3K4me3"
jmark <- "H3K27me3"
inf.mat <- file.path(hubprefix, paste0("jyeung/data/scChiC/mara_analysis_BM-AllMerged2_Peaks_1000/", jmark, "/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, ".txt"))
outf.mat <- file.path(hubprefix, paste0("jyeung/data/scChiC/mara_analysis_BM-AllMerged2_Peaks_1000/", jmark, "/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, ".RemoveLargePeaks.txt"))
outf.mat.colscaled <- file.path(hubprefix, paste0("jyeung/data/scChiC/mara_analysis_BM-AllMerged2_Peaks_1000/", jmark, "/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, ".RemoveLargePeaks.ColumnScaled.txt"))
outf.mat.rowscaled <- file.path(hubprefix, paste0("jyeung/data/scChiC/mara_analysis_BM-AllMerged2_Peaks_1000/", jmark, "/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, ".RemoveLargePeaks.RowScaled.txt"))
outpdf <- file.path(hubprefix, paste0("jyeung/data/scChiC/mara_analysis_BM-AllMerged2_Peaks_1000/", jmark, "/mara_input/sitecount_mats/hiddenDomains_motevo_merged.closest.long.scale_0.center_0.byrow_0.", jmark, ".RemoveLargePeaks.pdf"))

# assertthat::assert_that(!file.exists(outf.mat))
# assertthat::assert_that(!file.exists(outf.mat.colscaled))
# assertthat::assert_that(!file.exists(outf.mat.rowscaled))
# assertthat::assert_that(!file.exists(outpdf))


mat <- read.table.handlerows(inf.mat)

inf.exprs <- file.path(hubprefix, paste0("jyeung/data/scChiC/mara_analysis_BM-AllMerged2_Peaks_1000/", jmark, "/mara_input/count_mats_peaks_from_GLMPCA/BM_", jmark, ".BM_AllMerged.VarCorrection.Poisson.GLMPCA.txt"))
outf.exprs <- file.path(hubprefix, paste0("jyeung/data/scChiC/mara_analysis_BM-AllMerged2_Peaks_1000/", jmark, "/mara_input/count_mats_peaks_from_GLMPCA/BM_", jmark, ".BM_AllMerged.VarCorrection.Poisson.GLMPCA.RemoveLargePeaks.txt"))
# assertthat::assert_that(!file.exists(outf.exprs))


E <- read.table.handlerows(inf.exprs)

# Get common rows  --------------------------------------------------------

rows.common <- intersect(rownames(mat), rownames(E))


# Get distances -----------------------------------------------------------

coords.dat <- data.frame(coord = rows.common,
                         chromo = sapply(rows.common, GetChromo), 
                         jstart = as.numeric(sapply(rows.common, GetStart)), 
                         jend = as.numeric(sapply(rows.common, GetEnd)),
                         stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(jdist = (jend - jstart),
         jdist.log10 = log10(jdist))

pdf(outpdf, useDingbats = FALSE)

ggplot(coords.dat, aes(x = jdist)) + 
  geom_density() 

ggplot(coords.dat, aes(x = log10(jdist))) + 
  geom_density() +
  geom_vline(xintercept = 4.6)


# Filter distance ---------------------------------------------------------

jthres <- 0.5
mm <- mixtools::normalmixEM(x = coords.dat$jdist.log10)

xline <- min(mm$x[which(mm$posterior[, 1] < jthres)])
print(paste("Predicted xline:", xline))

# if (xline < 4){
  # print(paste("xline less than 4, setting xline to 4..."))
  print(paste("setting xline to 4..."))
  xline <- 4
# }

ggplot(coords.dat, aes(x = log10(jdist))) + 
  geom_density() +
  geom_vline(xintercept = xline)

plot(density(coords.dat$jdist.log10), main = jthres, xlab = "log10 Distance")
abline(v = xline, col = 'blue')
plot.mixEM(mm, whichplots = 2, xlab2 = "log10 Distance", main2 = paste(signif(xline, digits = 2), "thres:", jthres))
abline(v = xline, col = 'blue')

dev.off()

# Get clean mat  ----------------------------------------------------------

rows.keep <- subset(coords.dat, jdist.log10 <= xline)$coord

length(rows.common)
length(rows.keep)



# Write new sitecount mat and expr mat ------------------------------------

Nmat <- mat[rows.keep, ]
N <- data.frame(Gene.ID = rows.keep, Nmat, stringsAsFactors = FALSE)
N.colscaled <- data.frame(Gene.ID = rownames(Nmat), as.matrix(scale(Nmat, center = FALSE, scale = TRUE)), stringsAsFactors = FALSE)
N.rowscaled <- data.frame(Gene.ID = rownames(Nmat), as.matrix(t(scale(t(Nmat), center = FALSE, scale = TRUE))), stringsAsFactors = FALSE)

Enew <- data.frame(Gene.ID = rows.keep, E[rows.keep, ], stringsAsFactors = FALSE)

# Write new  --------------------------------------------------------------

fwrite(N, file = outf.mat, sep = "\t")
fwrite(N.colscaled, file = outf.mat.colscaled, sep = "\t")
fwrite(N.rowscaled, file = outf.mat.rowscaled, sep = "\t")
fwrite(Enew, file = outf.exprs, sep = "\t")



