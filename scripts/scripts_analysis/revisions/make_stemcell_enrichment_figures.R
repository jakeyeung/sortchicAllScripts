# Jake Yeung
# Date of Creation: 2019-11-17
# File: ~/projects/scchic/scripts/scripts_analysis/revisions/make_stemcell_enrichment_figures.R
# Show all the data that we got 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(scchicFuncs)


# function ---------------------------------------------------

LoadDatsCalculateVars <- function(out.lda, out.lda.predict, jsettings){
  # do stuff here
  
  umap.out <- umap(posterior(out.lda)$topics, config = jsettings)
  dat.umap.long.lda <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2])
  
  dat.umap.long.lda$is.stem <- FALSE
  
  dat.umap.long.lda <- DoLouvain(topics.mat = posterior(out.lda)$topics, custom.settings.louv = jsettings, dat.umap.long = dat.umap.long.lda)
  
  dat.umap.long.lda$louvain <- paste(dat.umap.long.lda$louvain, dat.umap.long.lda$is.stem, sep = "_")
  
  dat.pred <- predict(umap.out, out.lda.predict$topics)
  
  dat.pred.long <- data.frame(cell = rownames(dat.pred), umap1 = dat.pred[, 1], umap2 = dat.pred[, 2], stringsAsFactors = FALSE) %>%
    mutate(is.stem = TRUE)
  
  dat.pred.long <- DoLouvain(topics.mat = out.lda.predict$topics, custom.settings.louv = jsettings, dat.pred.long)
  dat.pred.long$louvain <- paste(dat.pred.long$louvain, dat.pred.long$is.stem, sep = "_")
  
  dat.umap.pred.merged <- bind_rows(dat.umap.long.lda, dat.pred.long)
  
  # calculate variance
  dat.impute.log <- log2(t(posterior(out.lda)$topics %*% posterior(out.lda)$terms) * jfac + jpseudo)
  dat.impute.log.enrich <- log2(t(out.lda.predict$topics %*% out.lda.predict$terms) * jfac + jpseudo)
  
  cells.var.chromo <- CalculateVarAll(dat.impute.log, jchromos)
  cells.var.chromo$experi <- "Control"
  cells.var.chromo.enrich <- CalculateVarAll(dat.impute.log.enrich, jchromos)
  cells.var.chromo.enrich$experi <- "StemCellEnriched"
  cells.var.chromo.merge <- rbind(cells.var.chromo, cells.var.chromo.enrich)
  
  umap.out.long.merge <- left_join(dat.umap.pred.merged, cells.var.chromo.merge)
  return(umap.out.long.merge)
}



# Settings ----------------------------------------------------------------

jsettings <- umap.defaults
jsettings$n_neighbors <- 25
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

jchromos <- c(paste("chr", seq(19), sep = ""), "chrX", "chrY")
jchromos.grep <- paste(jchromos, ":", sep = "")

jfac <- 10^6
jpseudo <- 0

outpdf <- paste0("/Users/yeung/data/scchic/pdfs/stemcell_analysis/debugging/projection_debugging.", Sys.Date(), ".pdf")

pdf(outpdf, useDingbats = FALSE)

# Load Zebrafish all  -----------------------------------------------------





# Init obj for bone marrow  -----------------------------------------------



dat.umap.all <- list()  # for the four marks 

# H3K4me1 analysis --------------------------------------------------------


# H3K4me1 is done in separate 
kchoose <- 50

inf.lda <- "/Users/yeung/data/scchic/from_cluster/oud3700/ldaAnalysisBins_PZ-Bl6-BM-StemCells/lda_out_meanfilt.PZ-Bl6-BM-StemCells_H3K4me1_matsMergedNonenriched_2019-09-29.CountThres0.K-30_35_50.Robj"
load(inf.lda, v=T)
kchoose.i <- which(sapply(out.lda, function(x) x@k) == kchoose)
out.lda <- out.lda[[kchoose.i]]

inf.proj <- "/Users/yeung/data/scchic/from_cluster/oud3700/ldaAnalysisBins_PZ-Bl6-BM-StemCells/projections/H3K4me1_matsMergedNonenriched.matsMergedEnriched.RData"
load(inf.proj, v=T)

umap.out.long.merge <- LoadDatsCalculateVars(out.lda, out.lda.predict, jsettings)

dat.umap.all$H3K4me1 <- umap.out.long.merge

ggplot(dat.umap.all$H3K4me1, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + 
  facet_wrap(~is.stem, ncol = 1) + 
  ggtitle("K4me1 by projection")

ggplot(dat.umap.all$H3K4me1, aes(x = umap1, y = umap2, color = is.stem)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_d(direction = -1) + 
  ggtitle("K4me1 by projection: variance")



# check all
inf.all <- "/Users/yeung/data/scchic/from_cluster/oud3700/ldaAnalysisBins_PZ-Bl6-BM-StemCells/lda_out_meanfilt.PZ-Bl6-BM-StemCells_H3K4me1_matsMergedAll_2019-09-29.CountThres0.K-30_35_50.OutObjs.RData"
load(inf.all, v=T)

jmat.k4me1 <- out.objs$tm.result$topics
umap.out <- umap(jmat.k4me1, config = jsettings)
dat.umap.long.merged.k4me1 <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2]) %>%
  rowwise() %>%
  mutate(is.stem = grepl("Bl6-BM-stem-cells", cell))

ggplot(dat.umap.long.merged.k4me1, aes(x = umap1, y = umap2, color = is.stem)) + 
  geom_point(alpha = 0.5) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("K4me1 stem cell by LDA all")

# plot boxplot 

dat.impute.log.k4me1 <- log2(t(out.objs$tm.result$topics %*% out.objs$tm.result$terms) * jfac + jpseudo)

cells.var.chromo.k4me1 <- CalculateVarAll(dat.impute.log.k4me1, jchromos) 

cells.var.chromo.k4me1 <- cells.var.chromo.k4me1 %>%
  rowwise() %>%
  mutate(is.stem = grepl("Bl6-BM-stem-cells", cell))

ggplot(cells.var.chromo.k4me1, aes(x = cell.var.within.sum.norm, fill = is.stem)) + geom_density(alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("K4me1 stem cell by LDA all")
  

dat.impute.log.k4me1.merge <- left_join(dat.umap.long.merged.k4me1, cells.var.chromo.k4me1)

ggplot(dat.impute.log.k4me1.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point(alpha = 0.5) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) + facet_wrap(~is.stem) + 
  ggtitle("K4me1 stem cell by LDA all")
  



# H3K4me3 -----------------------------------------------------------------

jmark <- "H3K4me3"

inf <- paste0("/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.", jmark, ".stringent_filter.RData")
assertthat::assert_that(file.exists(inf))
load(inf, v=T)

inf.proj <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-StemCells/PZ-Bl6-BM-StemCells_matsMerged_", jmark, "_2019-11-17.RData")
assertthat::assert_that(file.exists(inf.proj))
load(inf.proj, v=T)

umap.out.long.merge <- LoadDatsCalculateVars(out.objs$out.lda, out.lda.predict, jsettings)

dat.umap.all$H3K4me3 <- umap.out.long.merge

ggplot(dat.umap.all$H3K4me3, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_viridis_c(direction = -1) +
  facet_wrap(~is.stem, ncol = 1) + 
  ggtitle("K4me3 by projection")
# 
ggplot(dat.umap.all$H3K4me3, aes(x = umap1, y = umap2, color = is.stem)) + geom_point(alpha = 0.5) + theme_bw() +
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("K4me3 by projection")


# check when merged all

inf.merged <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-StemCells_part2/lda_outputs.nonenriched_enriched_merged2.binarize.TRUE.no_filt/lda_out_meanfilt.PZ-Bl6-BM-StemCells_matsMergedWithNonenriched_H3K4me3_2019-11-17.CountThres0.K-30_50.Robj"
load(inf.merged, v=T)

out.lda <- out.lda[[2]]
jmat <- posterior(out.lda)$topics
umap.out <- umap(jmat, config = jsettings)
dat.umap.long.merged <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2]) %>%
  rowwise() %>%
  mutate(is.stem = grepl("B6BMSC", cell))

ggplot(dat.umap.long.merged, aes(x = umap1, y = umap2, color = is.stem)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("K4me3 by LDA merged")

# check variances 
# calculate variance
dat.impute.log <- log2(t(posterior(out.lda)$topics %*% posterior(out.lda)$terms) * jfac + jpseudo)

cells.var.chromo <- CalculateVarAll(dat.impute.log, jchromos) 

cells.var.chromo <- cells.var.chromo %>%
  rowwise() %>%
  mutate(is.stem = grepl("B6BMSC", cell))
 
ggplot(cells.var.chromo, aes(x = cell.var.within.sum.norm, fill = is.stem)) + geom_density(alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("K4me3 by LDA merged")

dat.umap.long.merged.merge <- left_join(dat.umap.long.merged, cells.var.chromo)

ggplot(dat.umap.long.merged.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("K4me3 by LDA merged") + 
  scale_color_viridis_c(direction = -1)



# 
# # Check K27me3  -----------------------------------------------------------
# 
# inf.merged <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-StemCells_part2/lda_outputs.nonenriched_enriched_merged2.binarize.TRUE.no_filt/lda_out_meanfilt.PZ-Bl6-BM-StemCells_matsMergedWithNonenriched_H3K27me3_2019-11-17.CountThres0.K-30_50.Robj"
# load(inf.merged, v=T)
# 
# out.lda <- out.lda[[2]]
# jmat <- posterior(out.lda)$topics
# umap.out <- umap(jmat, config = jsettings)
# dat.umap.long.merged <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2]) %>%
#   rowwise() %>%
#   mutate(is.stem = grepl("B6BMSC", cell))
# 
# ggplot(dat.umap.long.merged, aes(x = umap1, y = umap2, color = is.stem)) + 
#   geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# # check variances 
# # calculate variance
# dat.impute.log <- log2(t(posterior(out.lda)$topics %*% posterior(out.lda)$terms) * jfac + jpseudo)
# 
# cells.var.chromo <- CalculateVarAll(dat.impute.log, jchromos) 
# 
# cells.var.chromo <- cells.var.chromo %>%
#   rowwise() %>%
#   mutate(is.stem = grepl("B6BMSC", cell))
# 
# ggplot(cells.var.chromo, aes(x = cell.var.within.sum.norm, fill = is.stem)) + geom_density(alpha = 0.5) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())




# Check linneg ------------------------------------------------------------

inf.linneg <- "/Users/yeung/data/scchic/from_cluster/lda_outputs_linneg/lda_out_meanfilt.PZ-Bl6-BM-Linneg_H3K4me3_binfilt_cellfilt.CountThres0.K-30_35_50_mindist_0.4_mindist_processed_lda.Rdata"
load(inf.linneg, v=T)

jmat.k4me3all <- out.objs$tm.result$topics
umap.out <- umap(jmat.k4me3all, config = jsettings)
dat.umap.long.merged.k4me3all <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2]) %>%
  rowwise() %>%
  mutate(is.stem = grepl("Linneg", cell))

ggplot(dat.umap.long.merged.k4me3all, aes(x = umap1, y = umap2, color = is.stem)) + 
  geom_point(alpha = 0.5) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("K4me3 lineage negative")

# plot boxplot 

dat.impute.log.k4me3all <- log2(t(out.objs$tm.result$topics %*% out.objs$tm.result$terms) * jfac + jpseudo)

cells.var.chromo.k4me3all <- CalculateVarAll(dat.impute.log.k4me3all, jchromos) 

cells.var.chromo.k4me3all <- cells.var.chromo.k4me3all %>%
  rowwise() %>%
  mutate(is.stem = grepl("Linneg", cell))

ggplot(cells.var.chromo.k4me3all, aes(x = cell.var.within.sum.norm, fill = is.stem)) + geom_density(alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("K4me3 lineage negative")

ggplot(cells.var.chromo.k4me3all, aes(x = cell.var.within.sum.norm, fill = is.stem)) + geom_density(alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle("K4me3 lineage negative")

# merge with umap 
dat.umap.long.merged.k4me3all.merge <- left_join(dat.umap.long.merged.k4me3all, cells.var.chromo.k4me3all)

ggplot(dat.umap.long.merged.k4me3all.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point(alpha = 0.5) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1) +
  facet_wrap(~is.stem) + 
  ggtitle("K4me3 lineage negative")




# Show zebrafish  ---------------------------------------------------------


jmark <- "H3K4me1"
inf.zf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_ZFbonemarrow/lda_outputs.ZF-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.binarize.TRUE.2019-11-04/lda_out_meanfilt.ZF-", jmark, "_pcutoff_0.CountThres0.K-30_35_50.Robj")
assertthat::assert_that(file.exists(inf.zf))

load(inf.zf, v=T)
out.lda.zf <- out.lda[[3]]

jmat.zf <- posterior(out.lda.zf)$topics
umap.out <- umap(jmat.zf, config = jsettings)

dat.umap.long.merged.zf <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2]) %>%
  rowwise() %>%
  mutate(is.stem = grepl("ZFWKMCD41plus", cell))

ggplot(dat.umap.long.merged.zf, aes(x = umap1, y = umap2, color = is.stem)) + 
  geom_point(alpha = 0.5) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jmark, "Zebrafish")

# check variances 
# calculate variance
dat.impute.log.zf <- log2(t(posterior(out.lda.zf)$topics %*% posterior(out.lda.zf)$terms) * jfac + jpseudo)

cells.var.chromo.zf <- CalculateVarAll(dat.impute.log.zf, jchromos) 

cells.var.chromo.zf <- cells.var.chromo.zf %>%
  rowwise() %>%
  mutate(is.stem = grepl("ZFWKMCD41plus", cell))

ggplot(cells.var.chromo.zf, aes(x = cell.var.within.sum.norm, fill = is.stem)) + geom_density(alpha = 0.5) + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())

dat.umap.long.merged.zf.merge <- left_join(dat.umap.long.merged.zf, cells.var.chromo.zf)

ggplot(dat.umap.long.merged.zf.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  geom_point(alpha = 0.5) + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(jmark, "Zebrafish") + 
  scale_color_viridis_c(direction = -1) + 
  facet_wrap(~is.stem)


dev.off()

# project identical cell?  

umap.out.check <- umap(out.lda$topics, config = jsettings)
plot(umap.out.check$layout[, 1], umap.out.check$layout[, 2])



# check LDA only stem cells -----------------------------------------------

inf.stem <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-StemCells_part2/lda_outputs.nonenriched_enriched_merged.binarize.TRUE.no_filt/lda_out_meanfilt.PZ-Bl6-BM-StemCells_matsMerged_H3K4me3_2019-11-17.CountThres0.K-30_50.Robj")
assertthat::assert_that(file.exists(inf.stem))
load(inf.stem, v=T)

umap.out <- umap(posterior(out.lda[[2]])$topics, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2])

ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() 


# Lineage neg of K4me1 ----------------------------------------------------

inf.k4me1.linneg <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-Linneg/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.TRUE.no_filt/lda_out_meanfilt.PZ-Bl6-BM-Linneg_H3K4me1_binfilt_cellfilt.CountThres0.K-30_35_50.Robj"
load(inf.k4me1.linneg, v=T)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(posterior(out.lda[[3]])$topics, config = jsettings)

dat.umap.long.k4me1.linneg <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2]) %>%
  rowwise() %>%
  mutate(is.stem = grepl("Linneg", cell))

ggplot(dat.umap.long.k4me1.linneg, aes(x = umap1, y = umap2, color = is.stem)) + geom_point() + 
  theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ggtitle("H3K4me1 lineage negative")
  






# umap.pred <- predict(umap.out, data = out.lda.predict$topics[rep(c(1,2), 100), ])
# plot(umap.pred[, 1], umap.pred[, 2])

# check the merged output

# inf.check <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-StemCells_part2/lda_outputs.meanfilt_999.cellmin_0.cellmax_999999.binarize.TRUE.no_filt/lda_out_meanfilt.PZ-Bl6-BM-StemCells_matsMerged_H3K4me3_2019-10-24.CountThres0.K-30_35_50.Robj"
# load(inf.check, v=T)
# kchoose.i <- which(sapply(out.lda, function(x) x@k) == kchoose)
# out.lda <- out.lda[[kchoose.i]]
# 
# umap.out <- umap(posterior(out.lda)$topics, config = jsettings)
# dat.umap.long.lda <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2])
# 
# dat.umap.long.lda$is.stem <- sapply(dat.umap.long.lda$cell, function(x) grepl())

# 
# 
# # H3k27me3 ----------------------------------------------------------------
# 
# 
# jmark <- "H3K27me3"
# 
# inf <- paste0("/Users/yeung/data/scchic/tables/bamlist_for_merging_build95_B6_stringent/dat_umap_long_with_louvain.", jmark, ".stringent_filter.RData")
# assertthat::assert_that(file.exists(inf))
# load(inf, v=T)
# 
# inf.proj <- paste0("/Users/yeung/data/scchic/from_cluster/lda_outputs_stemcells_from_projection/PZ-Bl6-BM-StemCells_matsMerged_", jmark, "_2019-10-24.RData")
# load(inf.proj, v=T)
# 
# umap.out.long.merge <- LoadDatsCalculateVars(out.objs$out.lda, out.lda.predict, jsettings)
# 


# H3K9me3? ----------------------------------------------------------------

# nothing yet



# Check zebrafish  --------------------------------------------------------




