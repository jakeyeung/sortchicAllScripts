# Jake Yeung
# Date of Creation: 2019-11-23
# File: ~/projects/scchic/scripts/scripts_analysis/stemcell_analysis/6-analyze_unenriched_linneg_stemcells_merged.R
# 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(hash)
library(igraph)
library(umap)

library(Matrix)
library(topicmodels)

library(scchicFuncs)

# Load LDA output ---------------------------------------------------------

jmark <- "H3K4me3"
jstr <- "Linneg"
jbin <- FALSE
inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All/lda_outputs.PZ-Bl6-BM-All_", jstr, ".", jmark, ".2019-11-22.K-50.binarize.", jbin, "/ldaOut.PZ-Bl6-BM-All_", jstr, ".", jmark, ".2019-11-22.K-50.Robj")

# inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All/lda_outputs.PZ-Bl6-BM-AllMerged.H3K4me3.2019-11-22.K-50.binarize.FALSE/ldaOut.PZ-Bl6-BM-AllMerged.H3K4me3.2019-11-22.K-50.Robj"
# inf <- "./lda_outputs.PZ-Bl6-BM-All_Stemcells.H3K4me3.2019-11-22.K-50.binarize.FALSE/ldaOut.PZ-Bl6-BM-All_Stemcells.H3K4me3.2019-11-22.K-50.Robj"

inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All/lda_outputs.PZ-Bl6-BM-All_Stemcells.H3K4me3.2019-11-22.K-50.binarize.FALSE/ldaOut.PZ-Bl6-BM-All_Stemcells.H3K4me3.2019-11-22.K-50.Robj"

inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All/lda_outputs.PZ-Bl6-BM-All_Linneg.H3K4me3.2019-11-22.K-50.binarize.FALSE/ldaOut.PZ-Bl6-BM-All_Linneg.H3K4me3.2019-11-22.K-50.Robj"

# inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All/lda_outputs.PZ-Bl6-BM-All_Unenriched.H3K4me3.2019-11-22.K-50.binarize.FALSE/ldaOut.PZ-Bl6-BM-All_Unenriched.H3K4me3.2019-11-22.K-50.Robj"

inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks/lda_outputs.PZ-Bl6-BM-All_Linneg.H3K4me1.2019-11-23.K-50.binarize.TRUE/ldaOut.PZ-Bl6-BM-All_Linneg.H3K4me1.2019-11-23.K-50.Robj"

inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks/lda_outputs.PZ-Bl6-BM-All_Linneg.H3K4me1.2019-11-23.K-50.binarize.TRUE/ldaOut.PZ-Bl6-BM-All_Linneg.H3K4me1.2019-11-23.K-50.Robj"

inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks/lda_outputs.PZ-Bl6-BM-All_Stemcells.H3K27me3.2019-11-23.K-50.binarize.FALSE/ldaOut.PZ-Bl6-BM-All_Stemcells.H3K27me3.2019-11-23.K-50.Robj"

inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks/lda_outputs.PZ-Bl6-BM-All_Unenriched.H3K4me1.2019-11-23.K-50.binarize.TRUE/ldaOut.PZ-Bl6-BM-All_Unenriched.H3K4me1.2019-11-23.K-50.Robj"

inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks/lda_outputs.PZ-Bl6-BM-All_UnenrichedXLinNeg.H3K4me1.2019-11-23.K-50.binarize.TRUE/ldaOut.PZ-Bl6-BM-All_UnenrichedXLinNeg.H3K4me1.2019-11-23.K-50.Robj"

inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks/lda_outputs.PZ-Bl6-BM-All_UnenrichedXStemCells.H3K4me1.2019-11-23.K-50.binarize.TRUE/ldaOut.PZ-Bl6-BM-All_UnenrichedXStemCells.H3K4me1.2019-11-23.K-50.Robj"

inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All/lda_outputs.PZ-Bl6-BM-AllMerged.H3K4me3.2019-11-22.K-50.binarize.FALSE/ldaOut.PZ-Bl6-BM-AllMerged.H3K4me3.2019-11-22.K-50.Robj"

inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-Linneg/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.FALSE.no_filt/lda_out_meanfilt.PZ-Bl6-BM-Linneg_H3K4me3_binfilt_cellfilt.CountThres0.K-30_35_50.Robj"

# inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-Linneg/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.FALSE.no_filt/lda_out_meanfilt.PZ-Bl6-BM-Linneg_H3K27me3_binfilt_cellfilt.CountThres0.K-30_35_50.Robj"
inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-Linneg/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.FALSE.no_filt/lda_out_meanfilt.PZ-Bl6-BM-Linneg_H3K9me3_binfilt_cellfilt.CountThres0.K-30_35_50.Robj"

inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-Linneg/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.TRUE.no_filt/lda_out_meanfilt.PZ-Bl6-BM-Linneg_H3K27me3_binfilt_cellfilt.CountThres0.K-30_35_50.Robj"

inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All/lda_outputs.PZ-Bl6-BM-All_Unenriched.H3K4me3.2019-11-22.K-50.binarize.FALSE/ldaOut.PZ-Bl6-BM-All_Unenriched.H3K4me3.2019-11-22.K-50.Robj"

inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All/lda_outputs.PZ-Bl6-BM-AllMerged.H3K4me3.2019-11-22.K-50.binarize.FALSE/ldaOut.PZ-Bl6-BM-AllMerged.H3K4me3.2019-11-22.K-50.Robj"
assertthat::assert_that(file.exists(inf))
load(inf, v=T)

if (length(out.lda) > 1){
  out.lda <- out.lda[[1]]
}
tm.result <- posterior(out.lda)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(tm.result$topics, config = jsettings)

dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "-"))

dat.umap.long <- dat.umap.long %>%
  mutate(experi = gsub("-G2", "", experi))
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#FF0000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.5, size = 3) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") + 
  scale_color_manual(values = cbPalette)


# do variance
jchromos <- c(paste("chr", seq(19), sep = ""), "chrX", "chrY")
jfac <- 10^6
jpseudo <- 0

dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms) * jfac + jpseudo)
cells.var.chromo <- CalculateVarAll(dat.impute.log, jchromos)
dat.umap.long.var <- left_join(dat.umap.long, cells.var.chromo, by = "cell")

ggplot(dat.umap.long.var, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + 
  facet_wrap(~experi)

dat.umap.long.var %>% 
  ggplot(., aes(x = cell.var.within.sum.norm, fill = experi)) + geom_density(alpha = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~experi, ncol = 1)

dat.umap.long.var %>% 
  ggplot(., aes(x = cell.var.within.sum.norm, fill = experi)) + geom_density(alpha = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

# dat.umap.long.var %>% 
#   ggplot(., aes(x = cell.var.within.sum.norm, fill = experi)) + geom_density(alpha = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_wrap(~experi, ncol = 1)

# visualize a region
# take entire chromosome 15
jchromo <- "^chr15"

bins.sub <- grepl(jchromo, rownames(dat.impute.log))
jsub <- dat.impute.log[bins.sub, ]

# do hi variance cell?
jcell <- (subset(dat.umap.long.var, experi == "B6-13W1-BM-H3K4me3") %>% arrange(desc(cell.var.within.sum.norm)))$cell[[1]]
(subset(dat.umap.long.var, experi == "B6-13W1-BM-H3K4me3") %>% arrange(desc(cell.var.within.sum.norm)))$cell.var.within.sum.norm[[1]]
plot(jsub[2500:3000, jcell], type = "l", main = "High Var")

jcell <- (subset(dat.umap.long.var, experi == "B6-13W1-BM-H3K4me3") %>% arrange(cell.var.within.sum.norm))$cell[[1]]
(subset(dat.umap.long.var, experi == "B6-13W1-BM-H3K4me3") %>% arrange(cell.var.within.sum.norm))$cell.var.within.sum.norm[[1]]
plot(jsub[2500:3000, jcell], type = "l", main = "Low Var")

jcell <- (subset(dat.umap.long.var, experi == "PZ-ChIC-B6BMSC-H3K4me3-G1") %>% arrange(desc(cell.var.within.sum.norm)))$cell[[1]]
(subset(dat.umap.long.var, experi == "B6-13W1-BM-H3K4me3") %>% arrange(desc(cell.var.within.sum.norm)))$cell.var.within.sum.norm[[1]]
plot(jsub[2500:3000, jcell], type = "l", main = "High Var")

jcell <- (subset(dat.umap.long.var, experi == "PZ-ChIC-B6BMSC-H3K4me3-G1") %>% arrange(cell.var.within.sum.norm))$cell[[1]]
(subset(dat.umap.long.var, experi == "B6-13W1-BM-H3K4me3") %>% arrange(cell.var.within.sum.norm))$cell.var.within.sum.norm[[1]]
plot(jsub[2500:3000, jcell], type = "l", main = "Low Var")

jcell <- (subset(dat.umap.long.var, experi == "PZ-Bl6-BM-Linneg-H3K4me3") %>% arrange(desc(cell.var.within.sum.norm)))$cell[[1]]
(subset(dat.umap.long.var, experi == "B6-13W1-BM-H3K4me3") %>% arrange(desc(cell.var.within.sum.norm)))$cell.var.within.sum.norm[[1]]
plot(jsub[2500:3000, jcell], type = "l", main = "High Var")

jcell <- (subset(dat.umap.long.var, experi == "PZ-Bl6-BM-Linneg-H3K4me3") %>% arrange(cell.var.within.sum.norm))$cell[[1]]
(subset(dat.umap.long.var, experi == "B6-13W1-BM-H3K4me3") %>% arrange(cell.var.within.sum.norm))$cell.var.within.sum.norm[[1]]
plot(jsub[2500:3000, jcell], type = "l", main = "Low Var")

# plot raw data
dim(count.mat.orig)

# all stem cells
jexperi <- "PZ-ChIC-B6BMSC-H3K4me3-"
jexperi <- "B6-13W1-BM-H3K4me3"

subset(dat.umap.long.var, )
cells.stem <- subset(dat.umap.long.var, startsWith(experi, jexperi))$cell
print(length(cells.stem))
  
raw.vec <- count.mat.orig[grepl(jchromo, rownames(count.mat.orig)), cells.stem] %>%
  rowSums()

plot(log2(raw.vec), type = "l")

# total variance across experis?
dat.umap.long.var %>% 
  ggplot(., aes(x = cell.var.within.sum.norm, fill = experi)) + geom_density(alpha = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  facet_wrap(~experi, ncol = 1)

dat.umap.long.var %>% 
  filter(experi %in% c("PZ-Bl6-BM-Linneg-H3K4me3", "B6-13W1-BM-H3K4me3")) %>%
  ggplot(., aes(x = cell.var.within.sum.norm, fill = experi)) + geom_density(alpha = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

