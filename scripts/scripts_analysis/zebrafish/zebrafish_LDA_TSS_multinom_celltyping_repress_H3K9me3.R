# Jake Yeung
# Date of Creation: 2019-11-15
# File: ~/projects/scchic/scripts/scripts_analysis/zebrafish/zebrafish_LDA_TSS_multinom_celltyping_repress_H3K9me3.R
# 

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

library(hash)
library(igraph)
library(umap)



# Constants ---------------------------------------------------------------

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
# zscore.cutoff <- 1
zscore.cutoff <- "top500"
jprefix <- "ZFWKM"
jprefix.proj <- "ZFWKMCD41plus"

jmark <- "H3K9me3"
winsize <- 100000L

winsizes <- c(50000L, 100000L)

# outdir <- "/Users/yeung/data/scchic/pdfs/zebrafish/celltyping"


# tf.vec <- c(TRUE, FALSE)
# jexps <- c(0.25, 0.5)

exponentiate.ref.mat <- FALSE
jnorm.vec <- TRUE
jexppower <- 0.5

Kstr <- "30"

# for (exponentiate.ref.mat in tf.vec){
#   for (jnorm.vec in tf.vec){
#     for (jexppower in jexps){

# jtitle <- paste(paste0("Mark_", jmark), paste0("winsize_", winsize), paste0("Zcutoff_", zscore.cutoff), paste0("Exp_", jexppower), paste0("expref_", exponentiate.ref.mat), paste0("normrefvec_", jnorm.vec), sep = ".")
# outpdf <- file.path(outdir, paste0(jtitle, ".pdf"))
# print(jtitle)

# Load data ---------------------------------------------------------------


inf.annot <- paste0("/Users/yeung/data/scchic/tables/gene_tss.winsize_", winsize, ".species_drerio.nochr.bed")

for (winsize in winsizes){
  # init
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisTSS_ZFbonemarrow.v2/lda_outputs.PZ-ChIC-", jprefix, "-", jmark, 
                ".winsize_", winsize, ".merged.K-", Kstr, ".binarize.FALSE/ldaOut.PZ-ChIC-", jprefix, "-", jmark, 
                ".winsize_", winsize, ".merged.K-", Kstr, ".Robj")
  assertthat::assert_that(file.exists(inf))
  load(inf, v=T)
  
  if (length(out.lda) > 1){
    out.lda <- out.lda[[1]]
  } 
  
  count.mat.long <- CollapseRowsByGene(count.mat, as.long=TRUE)
  
  
  # load tx data ------------------------------------------------------------
  
  # from make_tx_dataset_zebrafish_WKM.R
  inf.WKM <- "/Users/yeung/data/scchic/public_data/Zebrafish_WKM/Baron_et_al_pseudobulk_Zebrafish_WKM.rds"
  dat.bulk <- readRDS(inf.WKM)
  dat.bulk.mat <- dcast(subset(dat.bulk, select = c(gene, celltype, exprs)), gene ~ celltype, value.var = "exprs")
  rownames(dat.bulk.mat) <- dat.bulk.mat$gene; dat.bulk.mat$gene <- NULL
  
  # set up reference data
  
  if (is.character(zscore.cutoff)){
    dat.bulk.keep <- dat.bulk %>%
      group_by(celltype) %>%
      mutate(rnk = rank(-zscore)) %>%
      filter(rnk < 500)
  } else {
    dat.bulk.keep <- dat.bulk %>%
      group_by(gene) %>%
      filter(max(abs(zscore)) > zscore.cutoff)
  }
  
  genes.chic <- sapply(rownames(count.mat), function(x) strsplit(x, ";")[[1]][[2]])
  ref.genes.keep <- intersect(as.character(dat.bulk.keep$gene), genes.chic)
  
  print(length(ref.genes.keep))
  
  
  # do for main
  count.filt <- dcast(subset(count.mat.long, Var1 %in% ref.genes.keep), Var1 ~ Var2, value.var = "count")
  rownames(count.filt) <- count.filt$Var1; count.filt$Var1 <- NULL
  count.filt <- as.matrix(count.filt)[ref.genes.keep, ]  # keep same order!
  
  
  # Analyze LDA -------------------------------------------------------------
  
  # add gene name to the coordinates (got lost in mat to sparse mat pipeline)
  topics.mat <- posterior(out.lda)$topics
  terms.mat <- posterior(out.lda)$terms
  
  colnames(topics.mat) <- paste0("topic_", colnames(topics.mat))
  rownames(terms.mat) <- paste0("topic_", rownames(terms.mat))
  
  print(head(out.lda@terms))
  
  # colnames(terms.mat) <- sapply(colnames(terms.mat), function(x) annots.hash[[x]])
  
  jsettings <- umap.defaults
  jsettings$n_neighbors <- 30
  jsettings$min_dist <- 0.1
  jsettings$random_state <- 123
  
  jsettings.1d <- jsettings; jsettings.1d$n_components <- 1
  
  umap.out <- umap(topics.mat, config = jsettings)
  umap.out.1d <- umap(topics.mat, config = jsettings.1d)
  
  dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)
  dat.umap.long.1d <- data.frame(cell = rownames(umap.out$layout), umap1.1d = umap.out$layout[, 1], stringsAsFactors = FALSE)
  
  
  m.umap.blank <- ggplot(dat.umap.long, aes(x = umap1, y = umap2)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    ggtitle(paste("LDA with windows around TSS. Winsize:", winsize))
  
  # get variance??
  dat.impute.log <- log2(t(topics.mat %*% terms.mat))
  rownames(dat.impute.log) <- gsub(";", "_", rownames(dat.impute.log))
  
  # intrachromosomal variance doesnt make sense when doing TSS, do genome wide
  # jchromos <- paste("chr", seq(25), sep = "")
  jchromos <- c("")
  
  dat.var <- CalculateVarAll(dat.impute.log, jchromos)
  
  dat.umap.long.varmerge <- left_join(dat.umap.long, dat.var)
  dat.umap.long.varmerge <- left_join(dat.umap.long.varmerge, dat.umap.long.1d)
  
  PlotXYWithColor(dat.umap.long.varmerge, xvar = "umap1", yvar = "umap2", cname = "cell.var.within.sum.norm")
  
  ggplot(dat.umap.long.varmerge, aes(x = cell.var.within.sum.norm, y = umap1.1d)) + geom_point() 
  
  m.umap.var <- ggplot(dat.umap.long.varmerge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c()
  print(m.umap.var)
  
  
  # Find celltypes ----------------------------------------------------------
  
  # get topics
  topics.sum <- OrderTopicsByEntropy(posterior(out.lda), jquantile = 0.99) %>%
    rowwise() %>%
    mutate(topic = gsub("^X", "topic_", topic))
  
  # show a topic
  
  dat.umap.long.merge <- left_join(dat.umap.long, data.frame(cell = rownames(topics.mat), topics.mat, stringsAsFactors = FALSE))
  colnames(dat.umap.long.merge) <- gsub("^X", "topic_", colnames(dat.umap.long.merge))
  
  jtop <- topics.sum$topic[[1]]
  
  m.umap.top <- PlotXYWithColor(dat.umap.long.merge, xvar = "umap1", yvar = "umap2", cname = jtop)
  
  # what are the genes? 
  
  # Plot Genes Across Topic
  
  # make terms mat once
  terms.mat.long <- data.table::melt(terms.mat)
  colnames(terms.mat.long) <- c("topic", "term", "weight")
  terms.mat.long$term <- as.character(terms.mat.long$term)
  terms.mat.long <- terms.mat.long %>%
    rowwise() %>%
    mutate(gene = strsplit(term, ";")[[1]][[2]]) %>%
    group_by(topic) %>%
    mutate(rnk = rank(-weight)) %>%
    arrange(desc(weight))
  
  
  # make subset depending on the topic
  keeptop <- 50
  
  jsub.terms <- subset(terms.mat.long, topic == jtop & rnk <= keeptop)
  top.genes <- jsub.terms$gene
  
  m.exprs <- PlotPseudobulkZscore(subset(dat.bulk, gene %in% top.genes))
  m.top <- PlotDecreasingWeights(as.data.frame(jsub.terms), jtitle = paste("Topic:", jtop, "top", keeptop))
  
  print(m.umap.top)
  print(m.exprs)
  print(m.top)
  
  keeptop <- 200
  # for all topics
  outpdf <- paste0("/Users/yeung/data/scchic/pdfs/zebrafish/TSS_analysis/ZF_TSS.", jmark, ".winsize_", winsize, ".", Sys.Date(), ".pdf")
  pdf(outpdf, useDingbats = FALSE)
  
  print(m.umap.blank)
  print(m.umap.var)
  for (jtop in topics.sum$topic){
    print(jtop)
    
    m.umap.top <- PlotXYWithColor(dat.umap.long.merge, xvar = "umap1", yvar = "umap2", cname = jtop)
    jsub.terms <- subset(terms.mat.long, topic == jtop & rnk <= keeptop)
    top.genes <- jsub.terms$gene
    
    m.exprs <- PlotPseudobulkZscore(subset(dat.bulk, gene %in% top.genes))
    m.top <- PlotDecreasingWeights(as.data.frame(jsub.terms), jtitle = paste("Topic:", jtop, "top", keeptop))
    
    print(m.umap.top)
    print(m.exprs)
    print(m.top)
  }
  dev.off()
}


# 
# 
# # Set up and do fits ------------------------------------------------------
# 
# if (exponentiate.ref.mat){
#   probs.lst.filt <- SetUpProbs(2^dat.bulk.mat[ref.genes.keep, ], norm.vec = jnorm.vec)
# } else {
#   probs.lst.filt <- SetUpProbs(dat.bulk.mat[ref.genes.keep, ], norm.vec = jnorm.vec)
# }
# 
# # invert probs
# # probs.lst.filt <- lapply(probs.lst.filt, function(x) 1-x)
# 
# all.cells <- colnames(count.filt)
# 
# # run fits
# source("/Users/yeung/projects/scchicFuncs/R/MultinomFunctions.R")
# 
# LL.ctype.lst <- FitMultinomsFn(count.filt, all.cells, probs.lst.filt, fn = function(x) max(x) - x)
# 
# # Summarize fits ----------------------------------------------------------
# 
# LL.dat <- SummarizeMultinomFits(LL.ctype.lst, count.filt, all.cells) %>%
#   mutate(is.stem = FALSE)
# 
# LL.dat.all <- MultinomFitsToLong(all.cells, LL.ctype.lst, count.filt)
# 
# # add to umap 
# dat.umap.long.merge <- left_join(dat.umap.long, LL.dat, by = "cell")
# 
# m.celltypes <- ggplot(dat.umap.long.merge %>% filter(p.max > log(0)), aes(x = umap1, y = umap2, color = p.max)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   facet_wrap(~ctype.pred) + ggtitle(jtitle) + 
#   scale_color_viridis_c()
# print(m.celltypes)
# 
# 
# 
# 
# print(m.celltypes)
# 
# 
# 
# # Do Louvain  -------------------------------------------------------------
# 
# dat.umap.long.louv <- DoLouvain(topics.mat, custom.settings.louv = jsettings, dat.umap.long = dat.umap.long)
# dat.umap.long.louv.merge <- left_join(dat.umap.long.louv, LL.dat)
# dat.umap.long.louv.merge.all <- left_join(LL.dat.all, dat.umap.long.louv)
# 
# dat.umap.long.louv.proj <- DoLouvain(out.lda.predict$topics, custom.settings.louv = jsettings, dat.umap.long = dat.umap.long.proj)
# dat.umap.long.louv.merge.proj <- left_join(dat.umap.long.louv.proj, LL.dat.proj)
# dat.umap.long.louv.merge.all.proj <- left_join(LL.dat.all.proj, dat.umap.long.louv.proj)
# 
# dat.umap.long.louv.both <- DoLouvain(rbind(topics.mat, out.lda.predict$topics), custom.settings.louv = jsettings, dat.umap.long = dat.umap.long.both)
# dat.umap.long.louv.merge.both <- left_join(rbind(LL.dat %>% mutate(is.stem = FALSE), LL.dat.proj %>% mutate(is.stem = TRUE)), dat.umap.long.louv.both)
# dat.umap.long.louv.merge.all.both <- left_join(rbind(LL.dat.all %>% mutate(is.stem = FALSE), LL.dat.all.proj %>% mutate(is.stem = TRUE)), dat.umap.long.louv.both)
# 
# 
# 
# 
# 
# # m.both <- ggplot(dat.umap.long.louv.both, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + 
# #   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
# #   scale_color_manual(values = cbPalette)
# 
# # annotate celltypes for each louvain
# dat.louv.preds <- dat.umap.long.louv.merge.all %>%
#   group_by(louvain, ctype.pred) %>%
#   summarise(LL = sum(LL),
#             p = mean(p))
# 
# dat.louv.preds.proj <- dat.umap.long.louv.merge.all.proj %>%
#   group_by(louvain, ctype.pred) %>%
#   summarise(LL = sum(LL),
#             p = mean(p))
# 
# dat.louv.preds.both <- dat.umap.long.louv.merge.all.both %>%
#   group_by(louvain, ctype.pred) %>%
#   summarise(LL = sum(LL),
#             p = mean(p))
# 
# # louvain 5 is erythyrocytes
# 
# # predict celltypes by summing likelihooods
# ctype.pred.louv <- dat.louv.preds %>%
#   group_by(louvain) %>%
#   mutate(p.louv = SoftMax(LL)) %>%
#   dplyr::filter(p.louv == max(p.louv))
# ctype.pred.louv.proj <- dat.louv.preds.proj %>%
#   group_by(louvain) %>%
#   mutate(p.louv = SoftMax(LL)) %>%
#   dplyr::filter(p.louv == max(p.louv))
# ctype.pred.louv.both <- dat.louv.preds.both%>%
#   group_by(louvain) %>%
#   mutate(p.louv = SoftMax(LL)) %>%
#   dplyr::filter(p.louv == max(p.louv))
# 
# # reannotate louvains 
# dat.umap.long.louv.annot <- left_join(dat.umap.long.louv, ctype.pred.louv)
# dat.umap.long.louv.proj.annot <- left_join(dat.umap.long.louv.proj, ctype.pred.louv.proj)
# dat.umap.long.louv.both.annot <- left_join(dat.umap.long.louv.both, ctype.pred.louv.both)
# 
# 
# m.louv <- ggplot(dat.umap.long.louv.annot, aes(x = umap1, y = umap2, color = ctype.pred)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_manual(values = cbPalette)
# m.louv.proj <- ggplot(dat.umap.long.louv.proj.annot, aes(x = umap1, y = umap2, color = ctype.pred)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_manual(values = cbPalette)
# m.louv.both <- ggplot(dat.umap.long.louv.both.annot, aes(x = umap1, y = umap2, color = ctype.pred)) + geom_point() + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_manual(values = cbPalette) + facet_wrap(~is.stem)
# 
# 
# # check variance in the bins?
# 
# pdf(outpdf, useDingbats = FALSE)
# print(m.celltypes)
# print(m.celltypes.proj)
# print(m.celltypes.bar)
# print(m.louv)
# print(m.louv.proj)
# print(m.louv.both)
# dev.off()
# #     }
# #   }
# # }
# 
# 
# 
