# Jake Yeung
# Date of Creation: 2020-01-03
# File: ~/projects/scchic/scripts/scripts_analysis/debugging_rerun/check_LDA_hiddenDomains.R
# 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)

library(topicmodels)

library(scchicFuncs)

# get bins
inf.lda.all <- "/Users/yeung/data/scchic/robjs/B6_objs/LDA_objects_all_marks.Rdata"
load(inf.lda.all, v=T)

# get peaks 
# inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysis_fromHiddenDomains/ldaAnalysisHiddenDomains_1000_build95.withchr.cells_from_bin_analysis_2019-04-15/lda_outputs.meanfilt_10.cellmin_0.cellmax_9999999.binarize.FALSE/lda_out_meanfilt.PZ-BM-H3K27me3.CountThres0.K-25_50.Robj"
# inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysis_fromHiddenDomains/ldaAnalysisHiddenDomains_1000_build95.withchr.cells_from_bin_analysis_2019-04-15/lda_outputs.meanfilt_10.cellmin_0.cellmax_9999999.binarize.TRUE/lda_out_meanfilt.PZ-BM-H3K4me3.CountThres0.K-25_50.Robj"
# inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysis_fromHiddenDomains/ldaAnalysisHiddenDomains_1000_build95.withchr.cells_from_bin_analysis_2019-04-15/lda_outputs.meanfilt_10.cellmin_0.cellmax_9999999.binarize.TRUE/lda_out_meanfilt.PZ-BM-H3K4me3.CountThres0.K-25_50.Robj"
# inf <- "/hpc/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysis_CorrPeakFilt_1000_build95_B6_from_traj.cells_from_bin_analysis/lda_outputs.meanfilt_10.cellmin_0.cellmax_9999999.binarize.TRUE/lda_out_meanfilt.B6-H3K4me1.CountThres0.K-25.Robj"
inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysis_fromHiddenDomains/ldaAnalysis_CorrPeakFilt_1000_build95_B6_from_traj.cells_from_bin_analysis/lda_outputs.meanfilt_10.cellmin_0.cellmax_9999999.binarize.TRUE/lda_out_meanfilt.B6-H3K4me1.CountThres0.K-50.Robj"
load(inf, v=T)

out.lda.choose <- out.lda[[1]]
tm.result <- posterior(out.lda.choose)

# cnames.old <- rownames(tm.result$topics)
# cnames.new <- SwitchColnames(cnames.old, jsplit = "-")
rownames(tm.result$topics) <- unname(rownames(tm.result$topics))

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(tm.result$topics, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)

ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# add variance
dat.impute.log <- log2(t(out.objs$H3K4me1$tm.result$topics %*% out.objs$H3K4me1$tm.result$terms))
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var <- CalculateVarAll(dat.impute.log, jchromos)

# show variance on peaks umap 
dat.umap.long.merge <- left_join(dat.umap.long, dat.var)
dat.umap.long.merge$plate <- sapply(dat.umap.long.merge$cell, function(x) ClipLast(x, jsep = "_"))

ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point()  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)
ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point()  + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~plate)

library(Rtsne)
tsne.out <- Rtsne(X = tm.result$topics, pca = 1)
dat.tsne.long <- data.frame(cell = rownames(tm.result$topics), tsne1 = tsne.out$Y[, 1], tsne2 = tsne.out$Y[, 2])
dat.tsne.long <- left_join(dat.tsne.long, dat.var)

ggplot(dat.tsne.long, aes(x = tsne1, y = tsne2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggplot(dat.tsne.long, aes(x = tsne1, y = tsne2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c()



