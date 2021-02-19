# Jake Yeung
# Date of Creation: 2021-01-29
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K4me1-H3K9me3_reseq/16-load_raw_counts_filter_cells_for_LDA_k4_k9_dynamic_regions.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(scchicFuncs)
library(JFuncs)

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
hubprefix <- "/home/jyeung/hub_oudenaarden"

inf.rds <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/count_name.H3K4me1.k4_k9_dynamic_bins.2021-01-30.rds")
out.rds <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/mat_K4_regions_only/count_name.H3K4me1.k4_genes_only.2021-02-18.rds")
inf2.rds <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/count_name.H3K9me3.k4_k9_dynamic_bins.2021-01-30.rds")
out2.rds <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_tables.BM.k4_k9_dynamic_regions/mat_K4_regions_only/count_name.H3K9me3.k9_genes_only.2021-02-18.rds")
mat <- readRDS(inf.rds)
mat2 <- readRDS(inf2.rds)


# Filter for K4me1 regions only  ------------------------------------------

rnames.keep.i <- grepl(";NM", rownames(mat))
cnames.keep.i <- !duplicated(colnames(mat))
cnames.keep.i2 <- !duplicated(colnames(mat2))

mat.filt <- mat[rnames.keep.i, cnames.keep.i]
mat.filt2 <- mat2[!rnames.keep.i, cnames.keep.i2]

print(dim(mat))
print(dim(mat.filt))

print(dim(mat2))
print(dim(mat.filt2))


# remove empty cells
mat.filt <- mat.filt[, colSums(mat.filt) > 0]
mat.filt2 <- mat.filt2[, colSums(mat.filt2) > 0]

print(dim(mat.filt))
print(dim(mat.filt2))

saveRDS(object = mat.filt, file = out.rds)
saveRDS(object = mat.filt2, file = out2.rds)
