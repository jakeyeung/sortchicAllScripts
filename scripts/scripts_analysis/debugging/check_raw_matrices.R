# Jake Yeung
# Date of Creation: 2019-11-26
# File: ~/projects/scchic/scripts/scripts_analysis/debugging/check_raw_matrices.R
# Check why H3K4me3 is strange??


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)

library(scchicFuncs)
library(irlba)

# Read stuff --------------------------------------------------------------

# h3k4me3?
 

inf <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-All_allmarks/PZ-Bl6-BM-All_Unenriched.H3K4me1.2019-11-23.rds"

inf <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-All/PZ-Bl6-BM-All_Merged.H3K4me3.2019-11-22.rds"

inf <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-All_allmarks/PZ-Bl6-BM-All_Merged.H3K4me3.2019-11-23.rds"

inf <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-All_allmarks/PZ-Bl6-BM-All_Unenriched.H3K4me3.2019-11-23.rds"

inf.am <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-All_allmarks/PZ-Bl6-BM-All_Unenriched.H3K4me3.2019-11-23.rds"
inf.orig <- "/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-All/PZ-Bl6-BM-All_Unenriched.H3K4me3.2019-11-22.rds"

# infs <- 
dat1 <- readRDS(inf.am)
dat2 <- readRDS(inf.orig)


plot(rowSums(dat1), rowSums(dat2))


# 1 - nnzero(dat) / length(dat)

lsi.out <- RunLSI(as.matrix(dat))

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(lsi.out$u, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2])

ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle(inf)


