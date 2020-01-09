# Jake Yeung
# Date of Creation: 2019-12-09
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/explore_LDA_ZF.R
# Plot LDA for different marks ZF

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(hash)
library(igraph)
library(umap)

library(topicmodels)

library(scchicFuncs)

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(ggrepel)

# Constants ---------------------------------------------------------------


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3")

jbin <- "FALSE"
jsuff <- "AllMerged"
# outdir <- "/home/jyeung/data/from_rstudioserver/scchic/pdfs/"
# outdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/LDA_downstream_ZF"
# outdir <- "/home/jyeung/hpc/scChiC/from_rstudioserver/LDA_downstream_ZF"
outdir <- "/home/jyeung/data/from_rstudioserver/scchic/pdfs/LDA_downstream_ZF.2020-01-05"
dir.create(outdir)

jwin <- 100000L

zscore.cutoff <- "top_500"

# load Chloe data
jchromos <- paste("chr", seq(25), sep = "")
inf.annot <- paste0("/home/jyeung/hpc/databases/gene_tss/gene_tss.winsize_", jwin, ".species_drerio.bed")
assertthat::assert_that(file.exists(inf.annot))


normtype <- "_cpmnorm"
jdate <- "2019-12-10"

inf.WKM <- paste0("/home/jyeung/hpc/scChiC/public_data/Baron_et_al_2019_Zebrafish_WKM/from_macbook/Baron_et_al_pseudobulk_Zebrafish_WKM", normtype, ".", jdate, ".rds")
assertthat::assert_that(file.exists(inf.WKM))

# UMAP settings
jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

# load tx data ------------------------------------------------------------

# from make_tx_dataset_zebrafish_WKM.R
dat.bulk <- readRDS(inf.WKM)
dat.bulk$celltype <- as.factor(dat.bulk$celltype)

# dat.bulk.mat <- dcast(subset(dat.bulk, select = c(gene, celltype, exprs)), gene ~ celltype, value.var = "exprs")
# rownames(dat.bulk.mat) <- dat.bulk.mat$gene; dat.bulk.mat$gene <- NULL

# set up reference data

# if (is.character(zscore.cutoff)){
#   rnk.cutoff <- as.numeric(strsplit(zscore.cutoff, "_")[[1]][[2]])
#   assertthat::assert_that(!is.na(rnk.cutoff))
#   dat.bulk.keep <- dat.bulk %>%
#     group_by(celltype) %>%
#     mutate(rnk = rank(-zscore)) %>%
#     filter(rnk < rnk.cutoff)
# } else {
#   dat.bulk.keep <- dat.bulk %>%
#     group_by(gene) %>%
#     filter(max(abs(zscore)) > zscore.cutoff)
# }

# recalculate 
# dat.bulk.nothrombo <- subset(dat.bulk, celltype != "thrombocytes") %>%
#   group_by(gene) %>%
#   mutate(zscore = scale(exprs, center = TRUE, scale = TRUE))


# randomly sample genes
jgenes.rand <- sample(x = as.character(dat.bulk$gene), size = 150)

ggplot(subset(dat.bulk, gene %in% jgenes.rand), aes(x = celltype, y = zscore)) + geom_point() + geom_boxplot()



# Plot gene expression ----------------------------------------------------

# ggplot(dat.bulk %>% filter(gene == "hlfb"), aes(x = celltype , y = zscore)) + geom_point()
# ggplot(dat.bulk %>% filter(gene == "hlfa"), aes(x = celltype , y = zscore)) + geom_point()
# ggplot(dat.bulk %>% filter(gene == "hlfa"), aes(x = celltype , y = zscore)) + geom_point()
# ggplot(dat.bulk %>% filter(gene == "s100a10b"), aes(x = celltype , y = zscore)) + geom_point()
# ggplot(dat.bulk %>% filter(gene == "s100a10a"), aes(x = celltype , y = zscore)) + geom_point()

keepn <- 250

# keepnvec <- c(50, 200, 250, 500)
keepnvec <- c(100, 150)

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

# Load and plot  ----------------------------------------------------------

for (keepn in keepnvec){
  
  
  for (jmark in jmarks){
    
    print(jmark)
    outf <- paste0("ZF_LDA_output.", jmark, ".binarize_", jbin, ".", jsuff, ".", normtype, ".keepn_", keepn, ".final.pdf")
    outpdf <- file.path(outdir, outf)
    inf <- paste0("/home/jyeung/data/from_cluster/scchic/LDA_outputs_all/ldaAnalysisBins_ZFWKM_allmarks_mergedtagged_dedupfixed_redo/lda_outputs.ZF_", 
                  jsuff, "_", jmark, 
                  ".TAcutoff_0.5.countscutoff_500.binfilt_cellfilt.2019-12-05.K-30.binarize.", jbin, 
                  "/ldaOut.ZF_", jsuff, "_", jmark, ".TAcutoff_0.5.countscutoff_500.binfilt_cellfilt.2019-12-05.K-30.Robj")
    assertthat::assert_that(file.exists(inf))
    load(inf, v=T)
    
    tm.result <- posterior(out.lda)
    topics.mat <- tm.result$topics
    colnames(topics.mat) <- paste("topic", colnames(topics.mat), sep = "_")
    
    
    
    # Plot UMAP  --------------------------------------------------------------
    
    umap.out <- umap(topics.mat, config = jsettings)
    
    dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
      rowwise() %>%
      mutate(experi = ClipLast(cell),
             plate = ClipLast(cell, jsep = "_"))
    
    m.umap.split <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
      facet_wrap(~experi) + ggtitle(jmark)
    
    dat.umap.long <- DoLouvain(topics.mat = topics.mat, custom.settings.louv = jsettings, dat.umap.long = dat.umap.long)
    
    dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
    
    dat.var <- CalculateVarAll(dat.impute.log, jchromos)
    dat.cellsizes <- data.frame(cell = colnames(count.mat), cellsize = log10(colSums(count.mat) / 5), stringsAsFactors = FALSE)
    
    dat.merge <- left_join(dat.umap.long, dat.var)
    dat.merge <- left_join(dat.merge, dat.cellsizes)
    
    m.var <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
      scale_color_viridis_c(direction = -1) + ggtitle(jmark)
    m.var.split <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
      scale_color_viridis_c(direction = -1) + ggtitle(jmark) + facet_wrap(~experi)
    m.var.plate <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
      scale_color_viridis_c(direction = -1) + ggtitle(jmark) + facet_wrap(~plate)
    
    m.size <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cellsize)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
      scale_color_viridis_c(direction = 1) + ggtitle(jmark)
    m.size.experi <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cellsize)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
      scale_color_viridis_c(direction = 1) + ggtitle(jmark) + facet_wrap(~experi)
    m.size.plate <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = cellsize)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
      scale_color_viridis_c(direction = 1) + ggtitle(jmark) + facet_wrap(~plate)
    
    m.louv <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = louvain)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
      scale_color_manual(values = cbPalette) + ggtitle(jmark)
    m.louv.experi <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = louvain)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
      scale_color_manual(values = cbPalette) + ggtitle(jmark) + facet_wrap(~experi)
    m.louv.plate <- ggplot(dat.merge, aes(x = umap1, y = umap2, color = louvain)) + 
      geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())  + 
      scale_color_manual(values = cbPalette) + ggtitle(jmark) + facet_wrap(~plate)
      
    
    # Get entropy for each topics ---------------------------------------------
    
    # load bins
    terms.mat.tmp <- tm.result$terms
    # rownames(terms.mat.tmp) <- rownames(terms.mat.tmp)
    annot.out <- AnnotateBins2(terms.mat = tm.result$terms, top.thres = 0.995, inf.tss = inf.annot, txdb = TxDb.Drerio.UCSC.danRer11.refGene, annodb = "org.Dr.eg.db", chromos.keep = jchromos)
    annot.out$terms.annot <- annot.out$terms.annot %>%
      rowwise() %>%
      mutate(termgene = ifelse(is.na(termgene), "", termgene))
    annot.out$terms.annot$topic <- sapply(annot.out$terms.annot$topic, function(x) paste("topic_", x, sep = ""))
    annot.out$terms.annot$gene <- sapply(annot.out$terms.annot$termgene, function(x){
      jsplit <- strsplit(x, ";")[[1]]
      if (length(jsplit) > 1){
        return(jsplit[[2]])
      } else {
        return("")
      }
    })
    
    topics.sum <- OrderTopicsByEntropy(tm.result, jquantile = 0.99)
    topics.sum$topic <- gsub("^X", "topic_", topics.sum$topic)
    
    # show an island
    topics.mat.tmp <- as.data.frame(topics.mat)
    topics.mat.tmp$cell <- rownames(topics.mat)
    topics.mat.tmp <- left_join(dat.umap.long, topics.mat.tmp)
    
    pdf(outpdf, useDingbats = FALSE)
    print(m.umap.split)
    print(m.var)
    print(m.var.split)
    print(m.var.plate)
    print(m.size)
    print(m.size.experi)
    print(m.size.plate)
    print(m.louv)
    print(m.louv.experi)
    print(m.louv.plate)
    
    for (i in seq(nrow(topics.sum))){
      jtop <- topics.sum$topic[[i]]
      print(jtop)
      # plot UMAP
      m.umap <- PlotXYWithColor(topics.mat.tmp, xvar = "umap1", yvar = "umap2", cname = jtop, cont.color = TRUE, jtitle = jtop) + scale_color_viridis_c()
      
      jsub.terms <- subset(annot.out$terms.annot, topic == jtop) %>%
        dplyr::filter(rnk <= keepn)
      
      m.top <- jsub.terms %>%
        ungroup() %>%
        mutate(term = forcats::fct_reorder(term, dplyr::desc(weight))) %>%
        ggplot(aes(x = term, y = log10(weight), label = gene)) +
        geom_point(size = 0.25) +
        theme_bw(8) +
        geom_text_repel(size = 2.5, segment.size = 0.1, segment.alpha = 0.25) +
        theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.text.x = element_text(angle = 45, hjust = 1, size = 2)) +
        xlab("") + ylab("Log10 Bin Weight") +
        ggtitle(paste("Top peak weights for:", jtop))
      
      # check gene expression across genes
      gfilt <- unique(jsub.terms$gene)
      # m.ctype <- subset(dat.bulk, gene %in% gfilt) %>%
      m.ctype <- subset(dat.bulk, gene %in% gfilt) %>%
        # mutate(celltype = ) %>%
        ggplot(., aes(x = forcats::fct_reorder(celltype, dplyr::desc(zscore), .fun = median), y = zscore)) + geom_boxplot() + 
        geom_point() + 
        theme_bw(24) + 
        theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
        ggtitle(jtop)
      
      print(m.umap)
      print(m.top)
      print(m.ctype)
      
    }
    dev.off()
  }  
  
  
}

