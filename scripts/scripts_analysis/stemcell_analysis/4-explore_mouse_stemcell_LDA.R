# Jake Yeung
# Date of Creation: 2019-11-20
# File: ~/projects/scchic/scripts/scripts_analysis/stemcell_analysis/4-explore_mouse_stemcell_LDA.R
# Stem cell 

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
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)

library(forcats)
library(ggrepel)


# load stem cell  ---------------------------------------------------------

jmark <- "H3K27me3"
jbin <- FALSE
Kstr <- "K-50"
# inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-StemCells/PZ-Bl6-BM-StemCells_matsMerged_H3K27me3_2019-11-17.RData"
# inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-StemCells_part2/lda_outputs.nonenriched_enriched_merged2.binarize.FALSE.no_filt/lda_out_meanfilt.PZ-Bl6-BM-StemCells_matsMergedWithNonenriched_H3K27me3_2019-11-17.CountThres0.K-50.Robj"
# inf <- "/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-StemCells_part2/lda_outputs.nonenriched_enriched_merged2.binarize.TRUE.no_filt/lda_out_meanfilt.PZ-Bl6-BM-StemCells_matsMergedWithNonenriched_H3K4me3_2019-11-17.CountThres0.K-30_50.Robj"
inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_PZ-Bl6-BM-StemCells_part2/lda_outputs.nonenriched_enriched_merged2.binarize.", jbin, ".no_filt/lda_out_meanfilt.PZ-Bl6-BM-StemCells_matsMergedWithNonenriched_", jmark, "_2019-11-17.CountThres0.", Kstr, ".Robj")
load(inf, v=T)
inf.mat <- paste0("/Users/yeung/data/scchic/count_mat_binfilt_cellfilt_for_LDA_PZ-Bl6-BM-StemCells_part2/PZ-Bl6-BM-StemCells_matsMergedWithNonenriched_", jmark, "_2019-11-17.RData")
load(inf.mat, v=T)
out.objs <- LoadLDABins(jmark = jmark, jbin = FALSE, top.thres = 0.99, inf = inf, convert.chr20.21.to.X.Y = TRUE, choose.k = 50)


out.lda <- out.lda[[1]]
tm.result <- posterior(out.lda)

# what are the genes that show this batch effect??
# head(out.objs$regions.annotated)
terms.filt <- AnnotateTermsToNearestGene(tm.result)
terms.filt$topic <- paste("topic", terms.filt$topic, sep = "_")

# pretty cnames
colnames(tm.result$topics) <- paste("topic", colnames(tm.result$topics), sep = "_")

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
umap.out <- umap(tm.result$topics, config = jsettings)
dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2]) %>%
  rowwise() %>%
  mutate(is.stem = grepl("B6BMSC", cell))

ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = is.stem)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())


  
# variance argument??
# calculate variance

jchromos <- c(paste("chr", seq(19), sep = ""), "chrX", "chrY")
jfac <- 10^6
jpseudo <- 0

dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms) * jfac + jpseudo)

cells.var.chromo <- CalculateVarAll(dat.impute.log, jchromos)

dat.umap.long.var <- left_join(dat.umap.long, cells.var.chromo, by = "cell")

ggplot(dat.umap.long.var, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + 
  theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  scale_color_viridis_c(direction = -1)


# Find batch effects?   ---------------------------------------------------

topics.sum <- OrderTopicsByEntropy(tm.result, jquantile = 0.99) %>%
  mutate(topic = gsub("X", "topic_", topic))

dat.umap.long.merge <- left_join(dat.umap.long, data.frame(cell = rownames(tm.result$topics), tm.result$topics, stringsAsFactors = FALSE))

# jtop <- topics.sum$topic[[2]]
jtop <- "topic_27"  # high in non-enriched
jtop <- "topic_25"  # high in non-enriched
jtop <- "topic_19"  # high in non-enriched
jtop <- "topic_5"  # high in non-enriched
jtop <- "topic_29"  # high in enriched
jtop <- "topic_30"  # high in enriched


jtops <- topics.sum$topic
pdf(paste0("/Users/yeung/data/scchic/pdfs/stemcell_analysis/debugging/topics_", jmark, ".pdf"), useDingbats = FALSE)
for (jtop in jtops){
  print(jtop)
  m <- PlotXYWithColor(dat.umap.long.merge, xvar = "umap1", yvar = "umap2", cname = jtop, jtitle = jtop)
  print(m)
  m.genes <- ggplot(subset(terms.filt, topic == jtop & rnk <= 100), aes(x = forcats::fct_reorder(term, weight, .desc = TRUE), y = weight, label = gene)) + 
    geom_point(size = 0.5) + 
    theme_bw(6) + 
    theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +  
    geom_text_repel(segment.alpha = 0.2) + 
    xlab("") + ggtitle(jtop)
  print(m.genes)
}
dev.off()


# 
# 
# jtop <- "topic_26"
# 
# 
# jtop <- "topic_26"
# jtop <- "topic_24"
# jtop <- "topic_27"
# subset(terms.filt, topic == jtop)
# jterm <- "chr11:51820000-51920000"  # Jade2
# jterm <- "chr7:115960000-116060000"  # Sox6
# jterm <- "chr11:98140000-98240000"  # "Med1"
# 
# # How many reads are we talking, in terms of raw reads?  ------------------
# 
# count.long <- melt(as.matrix(count.dat$counts))
# count.long$is.stem <- grepl("B6BMSC", count.long$Var2)
# 
# jsub <- subset(count.long, Var1 == jterm) %>%
#   dplyr::rename(term = Var1, cell = Var2, count = value) %>%
#   left_join(., dat.umap.long)
# # add umap coordinates
# 
# ggplot(jsub, aes(x = is.stem, y = count)) + geom_boxplot(outlier.size = NA, outlier.shape = NA) +
#   geom_jitter(width = 0.2) + 
#   ggtitle(jterm)
# 
# ggplot(jsub, aes(x = umap2, y = count)) + 
#   geom_jitter(height = 0.2) + 
#   ggtitle(jterm)
# 
# 
# jmod <- glm(value ~ is.stem, family = "binomial", data = jsub)
# 
# 
# # If we do celltyping... what do we get ?  --------------------------------
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
