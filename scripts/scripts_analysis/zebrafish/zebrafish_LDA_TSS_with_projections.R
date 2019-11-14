# Jake Yeung
# Date of Creation: 2019-11-13
# File: ~/projects/scchic/scripts/scripts_analysis/zebrafish/zebrafish_LDA_TSS_with_projections.R
# With projections


rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

library(hash)
library(igraph)
library(umap)

library(topicmodels)

library(scchicFuncs)

library(ggrepel)


# Load data ---------------------------------------------------------------

# jmark <- "H3K4me1"
# # winsize <- 100000L
# winsize <- 50000L

jprefix <- "ZFWKM"
jprefix.proj <- "ZFWKMCD41plus"

jmark <- "H3K9me3"
winsize <- 50000L


inf.annot <- paste0("/Users/yeung/data/scchic/tables/gene_tss.winsize_", winsize, ".species_drerio.nochr.bed")

# init
inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisTSS_ZFbonemarrow.v2/lda_outputs.PZ-ChIC-", jprefix, "-", jmark, 
              ".winsize_", winsize, ".merged.K-30.binarize.FALSE/ldaOut.PZ-ChIC-", jprefix, "-", jmark, 
              ".winsize_", winsize, ".merged.K-30.Robj")
inf.proj <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisTSS_ZFbonemarrow.v2/lda_outputs.PZ-ChIC-", jprefix, "-", jmark, 
                   ".winsize_", winsize, ".merged.K-30.binarize.FALSE/projections/ldaOut.PZ-ChIC-", jprefix, "-", jmark, ".winsize_", winsize, ".merged.K-30.x.PZ-ChIC-", jprefix.proj, "-", jmark, ".winsize_", winsize, ".merged.RData")
assertthat::assert_that(file.exists(inf))
assertthat::assert_that(file.exists(inf.proj))

x <- load(inf, v=T)
if (length(out.lda) > 1){
  out.lda <- out.lda[[1]]
} 
y <- load(inf.proj, v=T)


# add gene name to the coordinates (got lost in mat to sparse mat pipeline)
topics.mat <- posterior(out.lda)$topics
terms.mat <- posterior(out.lda)$terms

colnames(topics.mat) <- paste0("topic_", colnames(topics.mat))
rownames(terms.mat) <- paste0("topic_", rownames(terms.mat))

print(head(out.lda@terms))

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jsettings.1d <- jsettings; jsettings.1d$n_components <- 1

umap.out <- umap(topics.mat, config = jsettings)
umap.out.1d <- umap(topics.mat, config = jsettings.1d)

dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2])
dat.umap.long.1d <- data.frame(cell = rownames(umap.out$layout), umap1.1d = umap.out$layout[, 1])

m.umap.blank <- ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste("LDA with windows around TSS. Winsize:", winsize))
print(m.umap.blank)

# get variance??
dat.impute.log <- log2(t(topics.mat %*% terms.mat))
rownames(dat.impute.log) <- gsub(";", "_", rownames(dat.impute.log))

dat.impute.log.proj <- log2(t(out.lda.predict$topics %*% out.lda.predict$terms))
rownames(dat.impute.log.proj) <- gsub(";", "_", rownames(dat.impute.log.proj))



# intrachromosomal variance doesnt make sense when doing TSS, do genome wide
jchromos <- c("")

dat.var <- CalculateVarAll(dat.impute.log, jchromos)
dat.var.proj <- CalculateVarAll(dat.impute.log.proj, jchromos)

dat.umap.long.varmerge <- left_join(dat.umap.long, dat.var)
dat.umap.long.varmerge <- left_join(dat.umap.long.varmerge, dat.umap.long.1d)

PlotXYWithColor(dat.umap.long.varmerge, xvar = "umap1", yvar = "umap2", cname = "cell.var.within.sum.norm")

ggplot(dat.umap.long.varmerge, aes(x = cell.var.within.sum.norm, y = umap1.1d)) + geom_point() 


print(m.umap.blank)

m.umap.var <- ggplot(dat.umap.long.varmerge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c()
print(m.umap.var)


# project on data
umap.out.proj <- predict(object = umap.out, data = out.lda.predict$topics)

dat.umap.long.proj <- data.frame(cell = rownames(umap.out.proj), umap1 = umap.out.proj[, 1], umap2 = umap.out.proj[, 2], stringsAsFactors = FALSE) %>%
  mutate(is.stem = TRUE)

dat.umap.long.merged <- rbind(dat.umap.long %>% mutate(is.stem = FALSE), dat.umap.long.proj)

ggplot(dat.umap.long.merged, aes(x = umap1, y = umap2, color = is.stem)) + geom_point()  + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


dat.var.merge <- rbind(dat.var %>% mutate(is.stem = FALSE), dat.var.proj %>% mutate(is.stem = TRUE))
dat.umap.long.varmerge.proj <- left_join(dat.umap.long.merged, dat.var.merge)


m.umap.var.proj <- ggplot(dat.umap.long.varmerge.proj, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + 
  theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c() + facet_wrap(~is.stem)
print(m.umap.var.proj)


