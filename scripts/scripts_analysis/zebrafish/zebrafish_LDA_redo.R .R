# Jake Yeung
# Date of Creation: 2019-11-24
# File: ~/projects/scchic/scripts/scripts_analysis/zebrafish/zebrafish_LDA_redo.R .R
# LDA redo 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)



# Load LDA  ---------------------------------------------------------------

# jmark <- "H3K4me3"
jmark <- "H3K4me1"
# jstr <- "Merged"
# jstr <- "Unenriched"
jstr <- "StemCells"
jbin <- "FALSE"

# inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZF_All/lda_outputs.PZ-ZF-All_Merged.", jmark, ".2019-11-22.K-50.binarize.TRUE/ldaOut.PZ-ZF-All_Merged.", jmark, ".2019-11-22.K-50.Robj")
inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZF_All/lda_outputs.PZ-ZF-All_", jstr, ".", jmark, ".2019-11-22.K-50.binarize.", jbin, "/ldaOut.PZ-ZF-All_", jstr, ".", jmark, ".2019-11-22.K-50.Robj")

# inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZF_All/lda_outputs.PZ-ZF-All_Merged.H3K4me3.2019-11-22.K-50.binarize.FALSE/ldaOut.PZ-ZF-All_Merged.H3K4me3.2019-11-22.K-50.Robj"
# inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZF_All/lda_outputs.PZ-ZF-All_Merged.H3K4me1.2019-11-22.K-50.binarize.TRUE/ldaOut.PZ-ZF-All_Merged.H3K4me1.2019-11-22.K-50.Robj"
# inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZF_All/lda_outputs.PZ-ZF-All_Merged.H3K27me3.2019-11-22.K-50.binarize.TRUE/ldaOut.PZ-ZF-All_Merged.H3K27me3.2019-11-22.K-50.Robj"
# inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZF_All/lda_outputs.PZ-ZF-All_Merged.H3K9me3.2019-11-22.K-50.binarize.TRUE/ldaOut.PZ-ZF-All_Merged.H3K9me3.2019-11-22.K-50.Robj"

# inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZFbonemarrow/lda_outputs.ZF-H3K4me1_pcutoff_0.CountThres0.K-30_35_50.binarize.FALSE/lda_out_meanfilt.ZF-H3K4me1_pcutoff_0.CountThres0.K-30_35_50.Robj"

assertthat::assert_that(file.exists(inf))

load(inf, v=T)

out.lda <- out.lda
tm.result <- posterior(out.lda)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(tm.result$topics, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2]) %>%
  rowwise() %>%
  mutate(is.stem = grepl("ZFWKMCD41plus", cell))

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = is.stem)) + geom_point(alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jmark)

# do some variance?

jchromos <- paste("chr", seq(25), sep = "")
jfac <- 10^6
jpseudo <- 0
dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms) * jfac + jpseudo)
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

dat.umap.long.merge <- left_join(dat.umap.long, dat.var)

ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point(alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jmark) + facet_wrap(~is.stem) + scale_color_viridis_c(direction = -1)

ggplot(dat.umap.long.merge, aes(x = cell.var.within.sum.norm, fill = is.stem)) + geom_density(alpha = 0.5) +
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 




