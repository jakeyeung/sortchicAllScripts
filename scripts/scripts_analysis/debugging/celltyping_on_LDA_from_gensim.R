# Jake Yeung
# Date of Creation: 2019-12-02
# File: ~/projects/scchic/scripts/scripts_analysis/debugging/check_LDA_dedupfixed.R
# LDA dedup fixed

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(hash)
library(igraph)
library(umap)
library(Seurat)
library(scchicFuncs)

library(preprocessCore)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(DESeq2)

library(ggrepel)

library(here)

library(topicmodels)

setwd(here())

# Load data  --------------------------------------------------------------

# jmark <- "H3K27me3"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

for (jmark in jmarks){
  
  # jmark <- "H3K27me3"
  # jsuff <- "UnenrichedXLinneg"
  # jsuff <- "Unenriched"
  # jsuff <- "Unenriched"
  # jsuff <- "Linneg"
  # jsuff <- "AllMerged"
  jsuff <- "UnenrichedXStemCells"
  jbin <- FALSE
  
  outpdf <- paste0("/Users/yeung/data/scchic/pdfs/stemcell_linneg_analysis.redo/", jmark, ".", jsuff, ".bin_", jbin, ".celltype_using_topics_from_gensim.pdf")
  
  if (file.exists(outpdf)){
    print(paste("Outpdf exists, skipping for", jmark))
    next
  }
  
  
  
  inf <- paste0("/Users/yeung/data/scchic/from_cluster/ldaAnalysisBins_B6BM_All_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.B6BM_", jsuff, "_", jmark, ".TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.binarize.", jbin, "/ldaOut.B6BM_", jsuff, "_", jmark, ".TAcutoff_0.5.countscutoff_1000.binfilt_cellfilt.2019-12-05.K-30.Robj")
  if (!file.exists(inf)){
    print(inf)
    print("Inf doesnt exist, skipping")
    next
  }
  assertthat::assert_that(file.exists(inf))
  
  load(inf, v=T)
  
  if (is.na(all(count.mat.orig))){
    count.mat.orig <- count.mat
  }
  
  # dat <- RunLSI(as.matrix(count.mat))
  # jsettings <- umap.defaults
  # jsettings$n_neighbors <- 30
  # jsettings$min_dist <- 0.1
  # jsettings$random_state <- 123
  # umap.out <- umap(dat$u, config = jsettings)
  # dat.umap.long.lsi <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)
  # ggplot(dat.umap.long.lsi, aes(x = umap1, y = umap2)) + geom_point() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
  tm.result <- posterior(out.lda)
  topics.mat <- tm.result$topics
  terms.mat <- tm.result$terms
  
  jsettings <- umap.defaults
  jsettings$n_neighbors <- 30
  jsettings$min_dist <- 0.1
  jsettings$random_state <- 123
  umap.out <- umap(topics.mat, config = jsettings)
  dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE)  %>%
    rowwise() 
  
  if (jmark == "H3K4me1"){
    dat.umap.long <- dat.umap.long %>%
      mutate(experi = strsplit(cell, "-")[[1]][[5]])  # needs adjusting
  } else {
    dat.umap.long <- dat.umap.long %>%
      mutate(experi = ClipLast(cell))  # needs adjusting
  }
  
  
  
  
  print(unique(dat.umap.long$experi))
  
  assertthat::assert_that(length(unique(dat.umap.long$experi)) < 6)  # need to double check experi
  
  ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.3) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
  
  ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.3) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~experi)
  
  dat.impute.log <- log2(t(topics.mat %*% terms.mat))
  jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
  dat.var <- CalculateVarAll(dat.impute.log, jchromos)
  
  
  # Do celltyping  ----------------------------------------------------------
  
  
  # Load bulk ---------------------------------------------------------------
  
  inf.bulkdat <- "/Users/yeung/data/scchic/public_data/E-MTAB-3079-query-results.fpkms.tsv"
  dat <- fread(inf.bulkdat, sep = "\t")
  
  colnames(dat) <- gsub(" ", "_", colnames(dat))
  dat.long <- gather(dat, key = "CellType", value = "FPKM", -c("Gene_ID", "Gene_Name")) %>%
    group_by(Gene_ID) %>%
    mutate(FPKM = replace_na(FPKM, 0)) %>%
    group_by(Gene_Name, CellType) %>%
    summarise(FPKM = sum(FPKM)) %>%
    rowwise() %>%
    mutate(logFPKM = log2(FPKM + 1))
  
  # normalize across samples?
  ggplot(dat.long, aes(x = CellType, y = logFPKM)) + geom_boxplot() +
    theme_bw() +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  dat.mat <- tidyr::spread(dat.long %>%
                             ungroup() %>%
                             # mutate(gene = paste(Gene_Name, Gene_ID, sep = ";")) %>%
                             mutate(gene = Gene_Name) %>%
                             dplyr::select(gene, CellType, logFPKM),
                           key = CellType, value = logFPKM)  %>%
    as.data.frame()
  rownames(dat.mat) <- dat.mat$gene; dat.mat$gene <- NULL
  
  cnames.tmp <- colnames(dat.mat)
  rnames.tmp <- rownames(dat.mat)
  dat.mat <- preprocessCore::normalize.quantiles(as.matrix(dat.mat), copy = TRUE)  # strong normalization,
  colnames(dat.mat) <- cnames.tmp
  rownames(dat.mat) <- rnames.tmp
  
  boxplot(dat.mat)
  
  dat.norm.long <- gather(data.frame(gene = rownames(dat.mat), dat.mat), key = "celltype", value = "exprs", -gene) %>%
    group_by(gene) %>%
    mutate(zscore = scale(exprs, center = TRUE, scale = TRUE))
  
  dat.norm.zscore.mat <- spread(dat.norm.long %>% dplyr::select(-exprs), key = "celltype", value = "zscore") %>%
    as.data.frame()
  rownames(dat.norm.zscore.mat) <- dat.norm.zscore.mat$gene
  dat.norm.zscore.mat$gene <- NULL
  
  annot.out <- AnnotateBins(terms.mat = terms.mat)
  
  terms.sum <- annot.out$terms.filt %>%
    group_by(gene) %>%
    dplyr::filter(rnk == min(rnk))
  
  term2gene <- hash(terms.sum$term, terms.sum$gene)
  gene2term <- hash(terms.sum$gene, terms.sum$term)
  
  genes.keep <- rownames(dat.mat)
  
  top.nterms <- 150
  terms.all <- unique(subset(terms.sum %>% arrange(desc(weight)), rnk <= top.nterms)$term)
  genes.all <- unlist(sapply(terms.all, function(x) term2gene[[x]]))
  
  
  # Do likelihoods ----------------------------------------------------------
  
  # create likelihoods
  jcutoff <- 1.8
  plot(density(rowMeans(dat.mat)))
  abline(v = jcutoff)
  
  dat.mat.filt <- dat.mat[apply(dat.mat, 1, max) > jcutoff, ]
  
  genes.keep <- intersect(genes.all, rownames(dat.mat.filt))
  terms.keep <- sapply(genes.keep, function(x) gene2term[[x]])
  
  dat.mat.filt <- dat.mat.filt[genes.keep, ]
  
  # redefine genes that we keep
  genes.filt <- intersect(rownames(dat.mat.filt), genes.keep)
  
  dat.mat.filt <- dat.mat.filt[genes.filt, ]
  
  probs.lst.filt <- SetUpProbs(dat.mat.filt, norm.vec = TRUE)
  
  # Prepare count mat -------------------------------------------------------
  
  count.filt <- count.mat.orig[terms.keep, ]
  rownames(count.filt) <- genes.keep
  
  
  # Run fits ----------------------------------------------------------------
  
  
  all.cells <- colnames(count.filt)
  names(all.cells) <- all.cells
  
  LL.ctype.lst <- FitMultinoms(count.filt, all.cells, probs.lst.filt, exppower = 0.25)
  
  
  # Summarize ---------------------------------------------------------------
  
  LL.dat <- SummarizeMultinomFits(LL.ctype.lst, count.filt, all.cells) %>%
    mutate(is.stem = grepl("Linneg", cell))
  
  m.celltypes.bar <- ggplot(LL.dat %>% group_by(is.stem, ctype.pred) %>% summarise(count = length(ctype.pred)) %>% group_by(is.stem) %>% mutate(total = sum(count)),
                            aes(x = ctype.pred, y = count / total, group = is.stem, fill = is.stem)) +
    geom_bar(position = "dodge", stat = "identity")  +
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  print(m.celltypes.bar)
  
  # Do UMAP  ----------------------------------------------------------------
  
  dat.merged <- left_join(dat.umap.long, LL.dat)
  dat.merged <- left_join(dat.merged, dat.var)
  
  dat.merged <- dat.merged %>%
    mutate(ctype.pred.stringent = ifelse(p.max > log(0.99), ctype.pred, NA))
  
  cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
  
  ggplot(dat.merged, aes(x = umap1, y = umap2, color = experi)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
    scale_color_manual(values = cbPalette)
  
  ggplot(dat.merged, aes(x = umap1, y = umap2, color = ctype.pred)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
    scale_color_manual(values = cbPalette)
  
  ggplot(dat.merged, aes(x = umap1, y = umap2, color = ctype.pred.stringent)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
    scale_color_manual(values = cbPalette, na.value="gray80") + facet_wrap(~ctype.pred)
  
  # show predicted celltypes have low variance
  ggplot(dat.merged, aes(x = cell.var.within.sum.norm)) + geom_density(alpha = 0.5, fill = 'blue')  + 
    theme_bw() + theme(aspect.ratio=0.1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
    scale_color_manual(values = cbPalette) + 
    facet_wrap(~ctype.pred, ncol = 1)
  
  # calculate average
  dat.var.avg <- subset(dat.merged, !is.na(ctype.pred.stringent)) %>%
    group_by(ctype.pred.stringent) %>%
    summarise(var.mean = mean(cell.var.within.sum.norm))
  
  ggplot(dat.var.avg, aes(x = forcats::fct_reorder(.f = ctype.pred.stringent, .x = var.mean, .fun = median, .desc = TRUE), y = var.mean)) + geom_bar(stat = "identity") + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  ggplot(dat.merged, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) + facet_wrap(~experi)
  
  ggplot(dat.merged, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) + facet_wrap(~ctype.pred.stringent)
  
  
  
  # Do celltyping based on topics  ------------------------------------------
  
  tm.result.tmp <- tm.result
  colnames(tm.result.tmp$topics) <- paste("Topic_", colnames(tm.result.tmp$topics), sep = "")
  
  topics.sum <- OrderTopicsByEntropy(tm.result.tmp, jquantile = 0.99)
  
  # add to dat
  dat.merged.topics <- left_join(dat.umap.long, data.frame(cell = rownames(tm.result.tmp$topics), tm.result.tmp$topics))
  
  terms.filt <- annot.out$terms.filt %>%
    mutate(topic = paste("Topic_", topic, sep = ""))
  
  # terms.filt <- data.frame(topic = paste0("Topic_", rownames(tm.result$terms)), as.data.frame(tm.result$terms)) %>%
  #   tidyr::gather(key = "term", value = "weight", -topic) %>%
  #   rowwise() %>%
  #   group_by(topic) %>%
  #   arrange(desc(weight)) %>%
  #   mutate(rnk = rank(-weight))
  # # add gene to terms name
  # terms.filt$gene <- sapply(terms.filt$term, function(x){
  #   xpretty <- paste0(strsplit(x, split = "\\.")[[1]][[1]], ":", strsplit(x, split = "\\.")[[1]][[2]], "-", strsplit(x, split = "\\.")[[1]][[3]])
  #   ifelse(!is.null(term2gene[[xpretty]]), term2gene[[xpretty]], xpretty)
  # })
  
  keeptop <- 150
  jtop <- "Topic_25"
  
  
  pdf(outpdf, useDingbats = FALSE)
  m.umap.experi <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.3) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  m.umap.experi2 <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point(alpha = 0.3) + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~experi)
  
  m.umap.var <- ggplot(dat.merged, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) + facet_wrap(~experi)
  m.umap.pred <- ggplot(dat.merged, aes(x = umap1, y = umap2, color = ctype.pred)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
    scale_color_manual(values = cbPalette)
  
  m.umap.pred2 <- ggplot(dat.merged, aes(x = umap1, y = umap2, color = ctype.pred.stringent)) + geom_point(alpha = 1, size = 5)  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right") +
    scale_color_manual(values = cbPalette)
  
  m.umap.stringent <- ggplot(dat.merged, aes(x = umap1, y = umap2, color = ctype.pred.stringent)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
    scale_color_manual(values = cbPalette, na.value="gray80") + facet_wrap(~ctype.pred)
  
  # show predicted celltypes have low variance
  m.umap.ctype.var <- ggplot(dat.merged, aes(x = cell.var.within.sum.norm)) + geom_density(alpha = 0.5, fill = 'blue')  + 
    theme_bw() + theme(aspect.ratio=0.1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom") +
    scale_color_manual(values = cbPalette) + 
    facet_wrap(~ctype.pred, ncol = 1)
  
  
  print(m.umap.experi)
  print(m.umap.var)
  print(m.umap.pred)
  print(m.umap.stringent)
  print(m.umap.ctype.var)
  
  for (jtop in topics.sum$topic){
    
    print(jtop)
    m.umap <- PlotXYWithColor(dat.merged.topics, xvar = "umap1", yvar = "umap2", cname = jtop)
    
    top.genes <- subset(terms.filt, topic == jtop & rnk <= keeptop)$gene
    assertthat::assert_that(length(top.genes) > 0)
    
    jsub <- subset(dat.norm.long, gene %in% top.genes)
    jsub.sorted.summarised <- jsub %>% group_by(celltype) %>% summarise(zscore = median(zscore)) %>% arrange(desc(zscore)) %>% dplyr::select(celltype)
    jlevels <- as.character(jsub.sorted.summarised$celltype)
    jsub$celltype <- factor(jsub$celltype, levels = jlevels)
    m.exprs <- ggplot(jsub,
                      aes(x = celltype , y = zscore)) +
      geom_boxplot(outlier.shape = NA) +
      # geom_violin() +
      geom_jitter(width = 0.1, size = 0.5) +
      # geom_line() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 4)) +
      ggtitle(paste(jtop, "Top:", keeptop, "N Unique Genes", length(top.genes)))
    # plot top 150 genes?
    jsub.terms <- subset(terms.filt, topic == jtop & rnk < keeptop) %>%
      ungroup() %>%
      mutate(term = forcats::fct_reorder(term, dplyr::desc(weight)))
    m.top <- jsub.terms %>%
      # mutate(term = forcats::fct_reorder(term, dplyr::desc(weight))) %>%
      ggplot(aes(x = term, y = log10(weight), label = gene)) +
      geom_point(size = 0.25) +
      theme_bw(8) +
      geom_text_repel(size = 2, segment.size = 0.1, segment.alpha = 0.25) +
      theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 3.5)) +
      xlab("") + ylab("Log10 Bin Weight") +
      ggtitle(paste("Top peak weights for:", jtop))
    
    # plot everything
    
    print(m.umap)
    print(m.exprs)
    print(m.top)
  }
  dev.off()
  
  
}

# 
# 
# # Plot raw data onto the UMAP  --------------------------------------------
# 
# # stem cell markers
# 
# dat.umap.louv <- DoLouvain(topics.mat = tm.result.tmp$topics, custom.settings.louv = jsettings, dat.umap.long = dat.umap.long %>% mutate(experi = gsub("-G2", "", experi)))
# 
# m.louv <- ggplot(dat.umap.louv, aes(x = umap1, y = umap2, color = louvain)) + geom_point(size = 5) + 
#   facet_wrap(~experi) + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
#   scale_color_manual(values = cbPalette)
# 
# # show smoothed topics 
# 
# jgenes.list <- c("Kit", "Hlf", "Meis1", "Etv6", "Ccl5", "Pid1", "Irf8", "Serpinb2", "Sirpa", "Pparg", "Mef2a", "Ikzf2", "Dhrs4", "Il6", "S100a8", "S100a9", "S100a7a", "Ebf1", "Ccl2", "Bach2", "Hbb-y", "F13a1")
# jterms.dat <- (subset(terms.filt, gene %in% jgenes.list) %>%
#   filter(weight == max(weight)))
# jgenes <- jterms.dat$gene
# jterms <- jterms.dat$term
# 
# jgenes.list[which(!jgenes.list %in% jgenes)]
# 
# jhash <- hash(jterms, jgenes)
# 
# # add to dat 
# mtmp <- 2^dat.impute.log[jterms, ]
# rownames(mtmp) <- sapply(rownames(mtmp), function(x) jhash[[x]])
# dtmp <- data.frame(cell = colnames(mtmp), t(mtmp))
# 
# # show raw data
# dat.umap.louv.exprsmerge <- left_join(dat.umap.louv, dtmp)
# 
# jgene <- "S100a7"
# jgene <- make.names(jgene)
# jterm <- grep(pattern = jgene, x = colnames(dat.umap.louv.exprsmerge), value = TRUE)
# PlotXYWithColor(dat.umap.louv.exprsmerge, xvar = "umap1", yvar = "umap2", cname = jterm)
# 
# 
# # Show enrichment in each louvain -----------------------------------------
# 
# dat.louv.sum <- dat.umap.louv %>%
#   group_by(louvain, experi) %>%
#   summarise(ncells = length(cell)) %>%
#   group_by(louvain) %>%
#   mutate(frac = ncells / sum(ncells))
# 
# m.bar <- ggplot(dat.louv.sum, aes(x = forcats::fct_reorder(.f = louvain, .x = frac, .fun = max, .desc = TRUE), y = frac, alpha = experi)) + geom_bar(stat = "identity", position = "dodge") + 
#   theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# 
# multiplot(m.louv, m.bar)

# 
# # Look at variance? -------------------------------------------------------
# 
# MergeCells <- function(mat, dat.split, jname, cname = "cell"){
#   cells.keep <- as.character(dat.split[[jname]][[cname]])
#   assertthat::assert_that(length(cells.keep) > 0)
#   return(rowSums(mat[, cells.keep]))
# }
# 
# print(m.louv)
# 
# # check the raw
# # for WT: see WT versus
# louv.highvar <- "3"
# louv.lowvar <- "2"
# jsplit <- split(dat.umap.louv, f = dat.umap.louv$louvain)
# 
# # # check granu traj 
# # jsplit.sub <- subset(dat.umap.louv, louvain == 5)
# # high.var <- as.character((subset(dat.var, cell %in% jsplit.sub$cell) %>% ungroup() %>% filter(cell.var.within.sum.norm > median(cell.var.within.sum.norm)))$cell)
# # low.var <- as.character((subset(dat.var, cell %in% jsplit.sub$cell) %>% ungroup() %>% filter(cell.var.within.sum.norm < median(cell.var.within.sum.norm)))$cell)
# 
# 
# x.prog.high.var <- MergeCells(count.mat.orig, jsplit, louv.highvar)
# x.prog.low.var <- MergeCells(count.mat.orig, jsplit, louv.lowvar)
# 
# # x.prog.high.var <- rowSums(count.mat.orig[, high.var])
# # x.prog.low.var <- rowSums(count.mat.orig[, low.var])
# 
# 
# # check chromo 15
# bins.filt <- grep("^chr15:", rownames(count.mat.orig), value = TRUE)
# # sort by numeri
# bins.filt.sorted <- bins.filt[order(sapply(bins.filt, GetStart), decreasing = FALSE)]
# 
# x1 <- x.prog.high.var[bins.filt.sorted]
# x2 <- x.prog.low.var[bins.filt.sorted]
# 
# x1 <- x1 / sum(x1)
# x2 <- x2 / sum(x2)
# 
# # x1 <- log2(x.prog.high.var[bins.filt.sorted] + 10)
# # x2 <- log2(x.prog.low.var[bins.filt.sorted] + 10)
# 
# # x1 <- log2(x1 / sum(x1) + 1)
# # x2 <- log2(x2 / sum(x2) + 1)
# 
# par(mfrow=c(2,1))
# plot(x1, type = "l", main = "High Var")
# plot(x2, type = "l", main = "Low Var")
# par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# 
# # do autocorrelation
# par(mfrow=c(1, 2))
# acf(x1, lag.max = 100, type = "correlation", demean = TRUE)
# acf(x2, lag.max = 100, type = "correlation", demean = TRUE)
# par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
# 
# 
# 
# 
# 
