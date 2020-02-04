# Jake Yeung
# Date of Creation: 2020-02-01
# File: ~/projects/scchic/scripts/rstudioserver_analysis/BM_all_merged/8-downstream_LDA_varfilt_GLMPCA.R
# Load LDA after filtering low var cells. Still do GLMPCA to correct for variance 

rm(list=ls())
jstart <- Sys.time()

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)

library(scchicFuncs)
library(hash)
library(igraph)
library(umap)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(ggrepel)


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
topn <- 150

# Load annotations --------------------------------------------------------




jinf.tss <- "/home/jyeung/hpc/databases/gene_tss/gene_tss_winsize.50000.bed"
assertthat::assert_that(file.exists(jinf.tss))


# inf.annot <- "/home/jyeung/hpc/scchic/public_data/Giladi_et_al_2018/giladi_pseudobulk_datsum_and_ncells.DESeq_and_QuantNorm.RData"
inf.annot <- "/home/jyeung/hpc/scChiC/public_data/Giladi_et_al_2018/giladi_pseudobulk_datsum_and_ncells.DESeq_and_QuantNorm.RData"
load(inf.annot, v=T)

dat.sum.long <- data.table::melt(dat.sum.norm.quantnorm)
colnames(dat.sum.long) <- c("gene", "celltype", "exprs")

dat.sum.long <- dat.sum.long %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE)) %>%
  filter(!is.nan(zscore))

jsuffix <- "KeepMorePlates"

infs <-  list("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-31.var_filt_keepPlates/lda_outputs.BM_H3K4me1.varcutoff_0.3.platesRemoved.K-30.binarize.FALSE/ldaOut.BM_H3K4me1.varcutoff_0.3.platesRemoved.K-30.Robj",
              "/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-31.var_filt_keepPlates/lda_outputs.BM_H3K4me3.varcutoff_0.3.KeepAllPlates.K-30.binarize.FALSE/ldaOut.BM_H3K4me3.varcutoff_0.3.KeepAllPlates.K-30.Robj",
              "/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-31.var_filt_keepPlates/lda_outputs.BM_H3K27me3.varcutoff_0.3.KeepAllPlates.K-30.binarize.FALSE/ldaOut.BM_H3K27me3.varcutoff_0.3.KeepAllPlates.K-30.Robj",
              "/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-31.var_filt_keepPlates/lda_outputs.BM_H3K9me3.varcutoff_0.3.KeepAllPlates.K-30.binarize.FALSE/ldaOut.BM_H3K9me3.varcutoff_0.3.KeepAllPlates.K-30.Robj")

names(infs) <- jmarks

for (jmark in jmarks){
  outname <- paste0("PZ_", jmark, ".topics_celltypes_Giladi.topn_", topn, ".", Sys.Date(), ".", jsuffix, ".pdf")
  outpdf <- file.path("/home/jyeung/hpc/scChiC/from_rstudioserver/pdfs_all/BM_LDA_downstream_topics_celltypes_Giladi", outname)
  print(jmark)
  print("Current time elapsed:")
  print(Sys.time() - jstart)
   
  # Load data  --------------------------------------------------------------
  
  # inf <- paste0("/home/jyeung/hpc/scChiC/raw_demultiplexed/LDA_outputs_all/ldaAnalysisBins_B6BM_All_allmarks.2020-01-31.var_filt/lda_outputs.BM_", jmark, ".varcutoff_0.3.platesRemoved.SmoothBinSize_1000.", jsuffix, ".K-30.binarize.FALSE/ldaOut.BM_", jmark, ".varcutoff_0.3.platesRemoved.SmoothBinSize_1000.", jsuffix, ".K-30.Robj")
  inf <- infs[[jmark]]
  load(inf, v=T)
  
  tm.result <- posterior(out.lda)
  colnames(tm.result$topics) <- paste("topic", colnames(tm.result$topics), sep = "_")
  rownames(tm.result$terms) <- paste("topic", rownames(tm.result$terms), sep = "_")
  
  
  topics.mat <- tm.result$topics
  
  
  # Plot it all -------------------------------------------------------------
  
  jsettings <- umap.defaults
  jsettings$n_neighbors <- 30
  jsettings$min_dist <- 0.1
  jsettings$random_state <- 123
  
  
  umap.out <- umap(topics.mat, config = jsettings)
  dat.umap.long <- data.frame(cell = rownames(umap.out[["layout"]]), umap1 = umap.out[["layout"]][, 1], umap2 = umap.out[["layout"]][, 2], stringsAsFactors = FALSE)
  dat.umap.long <- DoLouvain(topics.mat, jsettings, dat.umap.long)
  cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")
  ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = louvain)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_color_manual(values = cbPalette)
  
  
  # Plot variance -----------------------------------------------------------
  
  (jchromos <- paste("chr", c(seq(19)), sep = ""))
  dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
  
  dat.var <- CalculateVarAll(dat.impute.log, jchromos) %>%
    rowwise() %>%
    mutate(plate = ClipLast(as.character(cell), jsep = "_"))
  
  dat.var.merge <- left_join(dat.umap.long, dat.var)
  
  m.plates <- ggplot(dat.var.merge, aes(x = umap1, y = umap2, color = plate)) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_point() + scale_color_manual(values = cbPalette)
  
  m.var <- ggplot(dat.var.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_point() + scale_color_viridis_c(direction = -1)
  
  m.var.plates <- ggplot(dat.var.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    geom_point() + scale_color_viridis_c(direction = -1) + facet_wrap(~plate)
  
  # Do celltypes ------------------------------------------------------------
  
  
  
  topics.sum <- OrderTopicsByEntropy(tm.result, jquantile = 0.99)
  
  annots.out <- AnnotateBins2(tm.result$terms, top.thres = 0.995, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db")
  annots.out$terms.annot$gene <- sapply(annots.out$terms.annot$termgene, function(x) ifelse(is.na(x), NA, strsplit(x, ";")[[1]][[2]]))
  
  dat.topics <- data.frame(cell = rownames(topics.mat), topics.mat, stringsAsFactors = FALSE)
  
  dat.umap.long.merge <- left_join(dat.umap.long, dat.topics)
  
  # Plot a topikkukkc ------------------------------------------------------------
  
  jtopic <- topics.sum$topic[[1]]
  
  pdf(outpdf, useDingbats = FALSE)
  
  print(m.plates)
  print(m.var)
  print(m.var.plates)
  
  for (jtopic in topics.sum$topic){
    print(jtopic)
    m.umap <- PlotXYWithColor(dat.umap.long.merge, xvar = "umap1", yvar = "umap2", cname = jtopic)
    print(m.umap)
    terms.sub <- subset(annots.out$terms.annot, topic == jtopic)
    top.genes <- terms.sub$gene[1:topn]
    dat.sum.sub <- subset(dat.sum.long, gene %in% top.genes)
    m.exprs <- ggplot(dat.sum.sub,
                      aes(x = forcats::fct_reorder(celltype, zscore, .desc=TRUE), y = zscore)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.1, size = 0.5) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
      ggtitle(paste(jtopic, "Top:", topn, "N Unique Genes", length(top.genes))) + 
      xlab("")
    print(m.exprs)
    
    # plot top 150 genes?
    jsub.terms <- subset(terms.sub, topic == jtopic & rnk <= topn) %>%
      ungroup() %>%
      mutate(term = forcats::fct_reorder(term, dplyr::desc(weight)))
    m.top <- jsub.terms %>%
      ggplot(aes(x = term, y = log10(weight), label = gene)) +
      geom_point(size = 0.25) +
      theme_bw(8) +
      geom_text_repel(size = 2, segment.size = 0.1, segment.alpha = 0.25) +
      theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 3.5)) +
      xlab("") + ylab("Log10 Bin Weight") +
      ggtitle(paste("Top peak weights for:", jtopic))
    print(m.top)
  }
  dev.off()
}


print(Sys.time() - jstart)


