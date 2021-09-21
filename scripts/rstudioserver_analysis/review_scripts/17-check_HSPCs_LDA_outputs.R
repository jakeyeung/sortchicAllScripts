# Jake Yeung
# Date of Creation: 2021-09-07
# File: ~/projects/scchic/scripts/rstudioserver_analysis/review_scripts/17-check_HSPCs_LDA_outputs.R
# 

rm(list=ls())

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


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123
jsettings$spread <- 6

# Load LDA outputs --------------------------------------------------------

hubprefix <- "/home/jyeung/hub_oudenaarden"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

infs <- lapply(jmarks, function(jmark){
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_revisions/ldaAnalysis_HSPCsOnly_HSPCs_filt_countmat_bins_50kb.", jmark, ".2021-09-07.K-30/lda_outputs.HSPCs_filt_countmat_bins_50kb.", jmark, ".2021-09-07.K-30.binarize.FALSE/ldaOut.HSPCs_filt_countmat_bins_50kb.", jmark, ".2021-09-07.K-30.Robj"))
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

lda.outs <- lapply(infs, function(jinf){
  load(jinf, v=T)  # out.lda
  return(out.lda)
})

tm.lst <- lapply(lda.outs, function(jlda){
  tm <- posterior(jlda)
  tm <- AddTopicToTmResult(tm, jsep = "")
})



# Load metadata -----------------------------------------------------------

infs.meta <- lapply(jmarks, function(jmark){
  inf <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/cell_cluster_tables.batch_corrected_umaps.2020-12-28/rearranged_by_lineage/shuffled_cells/metadata_batch_corrected.arranged_by_lineage.shuffled.", jmark, ".2021-02-19.txt"))
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

dat.meta.lst <- lapply(infs.meta, fread)

# filter HSPCs
dat.meta.lst <- lapply(dats.meta, function(jdat){
  subset(jdat, cluster == "HSPCs", select = c(-umap1, -umap2, -louvain))
})



# Do umaps ----------------------------------------------------------------

dat.umap.lst <- lapply(tm.lst, function(tm){
  DoUmapAndLouvain(tm$topics, jsettings)
})

# show umap
m.lst <- lapply(jmarks, function(jmark){
  jdat <- dat.umap.lst[[jmark]]
  ggplot(jdat, aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    ggtitle(jmark) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

# annotate 



inf.annot <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Giladi_et_al_2018/giladi_pseudobulk_datsum_and_ncells.DESeq_and_QuantNorm.RData")
load(inf.annot, v=T)

dat.sum.long <- data.table::melt(dat.sum.norm.quantnorm)
colnames(dat.sum.long) <- c("gene", "celltype", "exprs")

dat.sum.long <- dat.sum.long %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE)) %>%
  filter(!is.nan(zscore))


jinf.tss <- file.path(hubprefix, "jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed")
assertthat::assert_that(file.exists(jinf.tss))
topn <- 150



# Check experimental and plate effects ------------------------------------

head(dat.umap.lst[[1]])
head(dat.umap.lst[[2]])




# check plate effects
dat.umap.annot.lst <- lapply(jmarks, function(jmark){
  dat.umap.lst <- dat.umap.lst[[jmark]] %>%
    left_join(., dat.meta.lst[[jmark]]) %>%
    rowwise() %>%
    mutate(experi = ClipLast(ClipLast(x = cell, jsep = "_"), jsep = "-"))
})

# check rep
cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")
m.batch <- lapply(jmarks, function(jmark){
  jdat <- dat.umap.annot.lst[[jmark]]
  ggplot(jdat, aes(x = umap1, y = umap2, color = batch)) + 
    geom_point() + 
    ggtitle(jmark) + 
    scale_color_manual(values = cbPalette) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})

m.rep <- lapply(jmarks, function(jmark){
  jdat <- dat.umap.annot.lst[[jmark]]
  ggplot(jdat, aes(x = umap1, y = umap2, color = jrep)) + 
    geom_point() + 
    ggtitle(jmark) + 
    scale_color_manual(values = cbPalette) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})


# check experi
m.experi <- lapply(jmarks, function(jmark){
  jdat <- dat.umap.annot.lst[[jmark]]
  ggplot(jdat, aes(x = umap1, y = umap2, color = experi)) + 
    geom_point() + 
    ggtitle(jmark) + 
    scale_color_manual(values = cbPalette) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
})





# Plot outputs ------------------------------------------------------------



for (jmark in jmarks){
  print(jmark)
  pdfout <- file.path(hubprefix, paste0("jyeung/data/scChiC/from_rstudioserver/post_submission/LDA_downstream_HSPCs/HSPCs_LDA.", jmark, ".", Sys.Date(), ".check_plates.pdf"))
  
  pdf(pdfout, useDingbats = FALSE)
  
  print(m.experi[[jmark]])
  
  tm.result <- tm.lst[[jmark]]
  
  
  topics.mat <- tm.result$topics
  topics.sum <- OrderTopicsByEntropy(tm.result, jquantile = 0.99)
  annots.out <- AnnotateBins2(tm.result$terms, top.thres = 0.995, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db")
  annots.out$terms.annot$gene <- sapply(annots.out$terms.annot$termgene, function(x) ifelse(is.na(x), NA, strsplit(x, ";")[[1]][[2]]))
  
  dat.topics <- data.frame(cell = rownames(topics.mat), topics.mat, stringsAsFactors = FALSE)
  
  dat.umap.long.merge <- dat.umap.lst[[jmark]] %>%
    left_join(., dat.topics, by = "cell")
  
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



