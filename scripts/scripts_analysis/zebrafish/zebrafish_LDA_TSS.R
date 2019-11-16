# Jake Yeung
# Date of Creation: 2019-11-12
# File: ~/projects/scchic/scripts/scripts_analysis/zebrafish/zebrafish_LDA_TSS.R
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
library(preprocessCore)
library(data.table)

library(here)

setwd(here())

GetGeneAnnotsHash <- function(inf.annot){
  dat.annot <- data.table::fread(inf.annot, col.names = c("chromo", "start", "end", "bname"))
  # add chr
  dat.annot$chromo <- paste("chr", dat.annot$chromo, sep = "")
  rnames.old <- paste(dat.annot$chromo, paste(dat.annot$start, dat.annot$end, sep = "-"), sep = ":")
  rnames.new <- dat.annot$bname
  annots.hash <- hash::hash(rnames.old, rnames.new)
}

AddGeneNameToRows <- function(mat, annots.hash){
  # mat rownmaes got stripped of gene names, add them back
  rnames.old <- rownames(mat)
  rnames.new <- sapply(rnames.old, function(x) annots.hash[[x]])
  rownames(mat) <- rnames.new
  return(mat)
}

PlotDecreasingWeights <- function(jsub.terms, jtitle = "", order.term.by.weight = TRUE, textsize = 3, themesize = 5){
  if (order.term.by.weight){
    jsub.terms <-  jsub.terms %>%
      mutate(term = forcats::fct_reorder(term, dplyr::desc(weight)))
  }
  m.top <- jsub.terms %>%
    ggplot(aes(x = term, y = log10(weight), label = gene)) +
    geom_point(size = 0.25) +
    theme_bw(themesize) +
    geom_text_repel(size = textsize, segment.size = 0.1, segment.alpha = 0.25) +
    theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 3.5)) +
    xlab("") + ylab("Log10 Bin Weight") +
    ggtitle(jtitle)
  return(m.top)
}

PlotPseudobulkZscore <- function(dat.bulk.sub, order.celltype.by.zscore = TRUE, xlabsize = 8, themesize = 12){
  if (order.celltype.by.zscore){
    dat.bulk.sub <- dat.bulk.sub %>% 
      mutate(celltype = forcats::fct_reorder(celltype, dplyr::desc(zscore), .fun = median))
  } 
  m.exprs <- ggplot(dat.bulk.sub,
                    aes(x = celltype , y = zscore)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 0.5) +
    theme_classic(themesize) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = xlabsize)) +
    ggtitle(paste("topic", jtop, "Top:", keeptop, "N Unique Genes", length(top.genes)))
  return(m.exprs)
  
}

# Load data ---------------------------------------------------------------

# jmark <- "H3K4me1"
# # winsize <- 100000L
# winsize <- 50000L

jprefix <- "ZFWKM"

jmarks <- c("H3K4me1", "H3K4me3", "H3K9me3")
winsizes <- c(50000L, 100000L)

norm.quants <- TRUE

# load tx data ------------------------------------------------------------

# from make_tx_dataset_zebrafish_WKM.R
inf.WKM <- "/Users/yeung/data/scchic/public_data/Zebrafish_WKM/Baron_et_al_pseudobulk_Zebrafish_WKM.rds"
dat.bulk <- readRDS(inf.WKM)


for (jmark in jmarks){
  for (winsize in winsizes){
    outpdf <- paste0("/Users/yeung/data/scchic/pdfs/zebrafish/TSS_analysis/ZF_TSS.", jmark, ".winsize_", winsize, ".", Sys.Date(), ".pdf")
    if (file.exists(outpdf)){
      print(paste(outpdf, "exists, skipping"))
      next
    }
    
    inf.annot <- paste0("/Users/yeung/data/scchic/tables/gene_tss.winsize_", winsize, ".species_drerio.nochr.bed")
    
    # init
    inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisTSS_ZFbonemarrow.v2/lda_outputs.PZ-ChIC-", jprefix, "-", jmark, 
                  ".winsize_", winsize, ".merged.K-30.binarize.FALSE/ldaOut.PZ-ChIC-", jprefix, "-", jmark, 
                  ".winsize_", winsize, ".merged.K-30.Robj")
    assertthat::assert_that(file.exists(inf))
    
    x <- load(inf, v=T)
    
    if (length(out.lda) > 1){
      out.lda <- out.lda[[1]]
    } 
    
    # add gene name to the coordinates (got lost in mat to sparse mat pipeline)
    topics.mat <- posterior(out.lda)$topics
    terms.mat <- posterior(out.lda)$terms
    
    colnames(topics.mat) <- paste0("topic_", colnames(topics.mat))
    rownames(terms.mat) <- paste0("topic_", rownames(terms.mat))
    
    print(head(out.lda@terms))
    
    # # if no gene names in rownames, then add them 
    # annots.hash <- GetGeneAnnotsHash(inf.annot)
    # count.mat <- AddGeneNameToRows(count.mat, annots.hash)
    
    # out.lda@terms <- sapply(out.lda@terms, function(x) annots.hash)
    
    
    # colnames(terms.mat) <- sapply(colnames(terms.mat), function(x) annots.hash[[x]])
    
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
    
    
    print(m.umap.blank)
    
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
}





