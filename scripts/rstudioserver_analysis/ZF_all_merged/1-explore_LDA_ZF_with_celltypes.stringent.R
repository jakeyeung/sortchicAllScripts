# Jake Yeung
# Date of Creation: 2020-04-11
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/1-explore_LDA_ZF_with_celltypes.R
# 

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

# UMAP settings
jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

ProcessLDA <- function(out.lda, jmark){
  print(jmark)
  tm.result <- posterior(out.lda)
  topics.mat <- tm.result$topics
  colnames(topics.mat) <- paste("topic", colnames(topics.mat), sep = "_")
  
  # Plot UMAP  --------------------------------------------------------------
  
  umap.out <- umap(topics.mat, config = jsettings)
  
  dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(experi = ClipLast(cell))
  
  m.umap.split <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~experi) + ggtitle(jmark)
  
  dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
  
  dat.var <- CalculateVarAll(dat.impute.log, jchromos)
  
  dat.merge <- left_join(dat.umap.long, dat.var)
  
  # Get entropy for each topics ---------------------------------------------
  
  # load bins
  terms.mat.tmp <- tm.result$terms
  # rownames(terms.mat.tmp) <- rownames(terms.mat.tmp)
  annot.out <- AnnotateBins2.ZF(terms.mat = tm.result$terms, top.thres = 0.995, inf.tss = inf.annot, txdb = TxDb.Drerio.UCSC.danRer11.refGene, annodb = "org.Dr.eg.db", chromos.keep = jchromos, skip.split = TRUE)
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
  outlst <- list()
  return(list(tm.result = tm.result, annot.pout = annot.out, topics.sum = topics.sum, dat.merge = dat.merge))
}

source("/home/jyeung/projects/scchic-functions/R/AuxLDA.R")

# Load LDA ----------------------------------------------------------------


jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks


# Constants ---------------------------------------------------------------

# keepn <- 250
keepnvec <- c(150)
cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

jwin <- 50000L
jprefix <- "/home/jyeung/hub_oudenaarden/jyeung"

# jsuffix <- ""
jsuffix <- ".imputevarfilt"
# jsuffix <- ".lessstringent"
indir <- file.path(jprefix, paste0("data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZF_AllMerged.winsize_", jwin, jsuffix))
outmain <- file.path(jprefix, paste0("data/zebrafish_scchic/from_rstudio/LDA_downstream"))
# outdir <- file.path(jprefix, paste0("data/scChiC/from_rstudioserver/LDA_downstream_ZF.", Sys.Date(), jsuffix))
outdir <- file.path(outmain, paste0("LDA_downstream_ZF.", Sys.Date(), jsuffix))
dir.create(outdir)

# outrdata <- file.path(outdir, paste0("dat_across_marks_outputs.RData"))

zscore.cutoff <- "top_500"

# load Chloe data
jchromos <- paste("chr", seq(25), sep = "")
inf.annot <- file.path(jprefix, paste0("data/databases/gene_tss/gene_tss.winsize_", jwin, ".species_drerio.bed"))
assertthat::assert_that(file.exists(inf.annot))


normtype <- "_cpmnorm"
normtype2 <- "_removethrombo"
jdate <- "2019-12-10"

inf.WKM <- file.path(jprefix, paste0("data/scChiC/public_data/Baron_et_al_2019_Zebrafish_WKM/from_macbook/Baron_et_al_pseudobulk_Zebrafish_WKM", normtype, ".", jdate, ".rds"))
assertthat::assert_that(file.exists(inf.WKM))

inf.WKMnothrombo <- file.path(jprefix, paste0("data/scChiC/public_data/Baron_et_al_2019_Zebrafish_WKM/from_macbook/Baron_et_al_pseudobulk_Zebrafish_WKM", normtype, ".", jdate, ".rds"))
assertthat::assert_that(file.exists(inf.WKMnothrombo))

# load tx data ------------------------------------------------------------

# from make_tx_dataset_zebrafish_WKM.R
dat.bulk <- readRDS(inf.WKM)
dat.bulk$celltype <- as.factor(dat.bulk$celltype)

dat.bulk.nothrombo <- readRDS(inf.WKMnothrombo)
dat.bulk.nothrombo$celltype <- as.factor(dat.bulk.nothrombo$celltype)


# Process -----------------------------------------------------------------


infs <- lapply(jmarks, function(jmark){
  if (jsuffix == ""){
    fname <- paste0("lda_outputs.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.Robj")
  } else if (jsuffix == ".imputevarfilt"){
    fname <- paste0("lda_outputs.counts_table_var_filt.", jmark, ".K-30.binarize.FALSE/ldaOut.counts_table_var_filt.", jmark, ".K-30.Robj")
  } else {
    # fname <- paste0("lda_outputs.counts_table_var_filt.", jmark, ".K-30.binarize.FALSE/ldaOut.counts_table_var_filt.", jmark, ".K-30.Robj")
    fname <- paste0("lda_outputs.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.Robj")
    # fname <- paste0()
  }
  inf <- file.path(indir, fname)
  print(inf)
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

print(infs)

out.lst <- lapply(infs, function(inf){
  load(inf, v=T)
  return(list(out.lda = out.lda, count.mat = count.mat))
})
 




# Plot gene expression ----------------------------------------------------




# Load and plot  ----------------------------------------------------------
 
 for (keepn in keepnvec){
   for (jmark in jmarks){
     print(jmark)
     
     out.lda <- out.lst[[jmark]]$out.lda
     count.mat <- out.lst[[jmark]]$count.mat
    
     outf <- paste0("ZF_LDA_output.", jmark, ".keepn_", keepn, ".final.pdf")
     outpdf <- file.path(outdir, outf)
     if (file.exists(outpdf)){
       print(paste(outpdf, "exists, skipping"))
       next
     }
     
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
     annot.out <- AnnotateBins2(terms.mat = tm.result$terms, top.thres = 0.995, inf.tss = inf.annot, txdb = TxDb.Drerio.UCSC.danRer11.refGene, annodb = "org.Dr.eg.db", chromos.keep = jchromos, skip.split = TRUE)
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
       
       dat.bulk.sub <- subset(dat.bulk, gene %in% gfilt)
       ngenes.kept <- length(unique(dat.bulk.sub$gene))
       m.ctype <- dat.bulk.sub %>%
         # mutate(celltype = ) %>%
         ggplot(., aes(x = forcats::fct_reorder(celltype, dplyr::desc(zscore), .fun = median), y = zscore)) + geom_boxplot() + 
         geom_jitter(width = 0.2) + 
         theme_bw(24) + 
         theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
         ggtitle(jtop, paste("Ngenes:", ngenes.kept))
       
       print(m.umap)
       print(m.top)
       print(m.ctype)
       
     }
     dev.off()
   }  
   
   
 }
 

