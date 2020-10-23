# Jake Yeung
# Date of Creation: 2020-08-17
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_with_others/integrate_BM_with_dbl_chic.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)


hubprefix <- "/home/jyeung/hub_oudenaarden"

# Load TSS ----------------------------------------------------------------

inf1.tss <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/Mouse_DblChIC_CTCF_K4me3_dChIC_run-1/tagged/count_tables_TSS/fBM-pCTCF-k4me3-i1.retagged.count_tables_TSS.csv")
inf2.tss <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/Mouse_DblChIC_CTCF_K4me3_dChIC_run-1/tagged/count_tables_TSS/fBM-pCTCF-k4me3-i2.retagged.count_tables_TSS.csv")

assertthat::assert_that(file.exists(inf1.tss))
assertthat::assert_that(file.exists(inf2.tss))

mat1.tss <- scchicFuncs::ReadMatTSSFormat(inf1.tss)
mat2.tss <- scchicFuncs::ReadMatTSSFormat(inf2.tss)

rowskeep <- intersect(rownames(mat1.tss), rownames(mat2.tss))
mats.merged <- cbind(mat1.tss[rowskeep, ], mat2.tss[rowskeep, ])
plot(density(colSums(mats.merged)))

cutmin <- 1000

mats.merged.filt <- mats.merged[, which(colSums(mats.merged) >= cutmin)]

# Load CTCF ---------------------------------------------------------------

inf1 <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/Mouse_DblChIC_CTCF_K4me3_dChIC_run-1/tagged/count_tables_CTCF/fBM-pCTCF-k4me3-i1.retagged.count_tables_CTCF.csv")
inf2 <- file.path(hubprefix, "jyeung/data/scChiC/raw_data_spikeins/Mouse_DblChIC_CTCF_K4me3_dChIC_run-1/tagged/count_tables_CTCF/fBM-pCTCF-k4me3-i2.retagged.count_tables_CTCF.csv")

assertthat::assert_that(file.exists(inf1))
assertthat::assert_that(file.exists(inf2))

mat1 <- scchicFuncs::ReadMatTSSFormat(inf1)
mat2 <- scchicFuncs::ReadMatTSSFormat(inf2)



# get CTCF counts per cell 
dat.ctcf.cuts1 <- data.frame(spikeincounts = colSums(mat1), cell = colnames(mat1), stringsAsFactors = FALSE)
dat.ctcf.cuts2 <- data.frame(spikeincounts = colSums(mat2), cell = colnames(mat2), stringsAsFactors = FALSE)

dat.ctcf.cuts <- as.data.frame(bind_rows(dat.ctcf.cuts1, dat.ctcf.cuts2))
rownames(dat.ctcf.cuts) <- dat.ctcf.cuts$cell

# Write mat to output -----------------------------------------------------

devmin <- 400
outf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/Mouse_DblChIC_CTCF_K4me3_dChIC_run-1/mats_for_LDA/count_mat_K36me3_CTCF_dbl_TSS.rds"
outf.filt <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/Mouse_DblChIC_CTCF_K4me3_dChIC_run-1/mats_for_LDA/count_mat_K36me3_CTCF_dbl_TSS.devmin_", devmin, ".rds")
outf.ctcf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_data_spikeins/Mouse_DblChIC_CTCF_K4me3_dChIC_run-1/mats_for_LDA/count_mat_K36me3_CTCF_dbl_CTCFcuts.rds"

saveRDS(object = mats.merged.filt, file = outf)
saveRDS(object = dat.ctcf.cuts, file = outf.ctcf)


# Filter rows  ------------------------------------------------------------

nvec <- dat.ctcf.cuts[colnames(mats.merged.filt), ]$spikeincounts
gdevs <- apply(mats.merged.filt, 1, function(xvec){
  scchicFuncs::binomial_deviance(x = xvec, p = sum(xvec) / sum(nvec), n = nvec)
})

plot(density(gdevs), main = "All genes", xlab = "Binomial Deviance")
abline(v = devmin)

genes.keep.vec <- gdevs[which(gdevs > devmin)]
genes.keep <- names(genes.keep.vec)

print(length(genes.keep))

mats.merged.filt2 <- mats.merged.filt[genes.keep, ]

saveRDS(object = mats.merged.filt2, file = outf.filt)


