# Jake Yeung
# Date of Creation: 2020-04-13
# File: ~/projects/scchic/scripts/rstudioserver_analysis/ZF_all_merged/3-GLMPCA_downstream.R
# Check GLMPCA outputs

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

library(JFuncs)

library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ChIPseeker)
library(GenomicRanges)


library(ggrepel)

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123


# cbPalette <- c("#696969", "#32CD32", "#F0E442", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

hubprefix <- "/home/jyeung/hub_oudenaarden"

keepn <- 150

jwin <- 50000L
jchromos <- paste("chr", seq(25), sep = "")

# Load pseudobulk data  ---------------------------------------------------

normtype <- "_cpmnorm"
jdate <- "2019-12-10"
inf.WKM <- file.path(hubprefix, paste0("jyeung/data/scChiC/public_data/Baron_et_al_2019_Zebrafish_WKM/from_macbook/Baron_et_al_pseudobulk_Zebrafish_WKM", normtype, ".", jdate, ".rds"))
assertthat::assert_that(file.exists(inf.WKM))

# from make_tx_dataset_zebrafish_WKM.R
dat.bulk <- readRDS(inf.WKM)
dat.bulk$celltype <- as.factor(dat.bulk$celltype)



# Iterate for all marks ---------------------------------------------------

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

for (jmark in jmarks){
  print(jmark)
  
  
  # Load GLMPCA output ------------------------------------------------------
  
  inmain <- paste0(hubprefix, "/jyeung/data/zebrafish_scchic/LDA_outputs/ldaAnalysisBins_ZF_AllMerged.winsize_50000")
  assertthat::assert_that(dir.exists(inmain))
  
  inf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/GLMPCA_outputs.ZF_AllMerged/ZF_", jmark, ".AllMerged.ZF_AllMerged.GLMPCA_var_correction.mergebinsize_1000.binskeep_500.covar_ncuts.var.log2.CenteredAndScaled.penalty_1.5.winsorize_TRUE.2020-04-13.RData")
  assertthat::assert_that(file.exists(inf))
  load(inf, v=T)
  
  
  # Load LDA ----------------------------------------------------------------
  
  infname <- paste0("lda_outputs.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.binarize.FALSE/ldaOut.count_mat.", jmark, ".countcutoffmin_1000-500-1000-1000.TAcutoff_0.5.chrfilt.bfilt.K-30.Robj")
  inf <- file.path(inmain, infname)
  assertthat::assert_that(file.exists(inf))
  
  load(inf, v=T)
  
  tm.result.lda <- posterior(out.lda)
  
  dat.umap.before <- DoUmapAndLouvain(tm.result.lda$topics, jsettings) %>%
    rowwise() %>%
    mutate(plate = ClipLast(cell, jsep = "_")) %>%
    dplyr::rename(louvain.before = louvain)
  
  # Show outputs ------------------------------------------------------------
  
  dat <- data.frame(cell = rownames(glm.out$factors), glm.out$factors, stringsAsFactors = FALSE)
  
  dat.umap <- DoUmapAndLouvain(glm.out$factors, jsettings) %>%
    rowwise() %>%
    mutate(plate = ClipLast(cell, jsep = "_")) %>%
    left_join(., subset(dat.umap.before, select = c(cell, louvain.before)))
  
  m.before <- ggplot(dat.umap.before, aes(x = umap1, y = umap2, color = louvain.before)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    scale_color_manual(values = cbPalette) 
  
  m.after.colorbybefore <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain.before)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    scale_color_manual(values = cbPalette) 
  
  m.after <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    scale_color_manual(values = cbPalette) 
  
  JFuncs::multiplot(m.before, m.after.colorbybefore, m.after, cols = 3)
  
  m.after.plates <- ggplot(dat.umap, aes(x = umap1, y = umap2, color = louvain)) + geom_point()  + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") + 
    scale_color_manual(values = cbPalette) + facet_wrap(~plate)
  print(m.after.plates)
  
  # Do celltyping ?  --------------------------------------------------------
  
  tm.result.glm <- list(topics = glm.out$factors, terms = t(glm.out$loadings))
  
  # Annotate terms ----------------------------------------------------------
  
  inf.annot <- file.path(hubprefix, paste0("jyeung/data/databases/gene_tss/gene_tss.winsize_", jwin, ".species_drerio.bed"))
  assertthat::assert_that(file.exists(inf.annot))
  
  annot.out <- AnnotateBins2(terms.mat = tm.result.glm$terms, top.thres = 0.995, inf.tss = inf.annot, txdb = TxDb.Drerio.UCSC.danRer11.refGene, annodb = "org.Dr.eg.db", chromos.keep = jchromos, skip.split = TRUE)
  annot.out$terms.annot <- annot.out$terms.annot %>%
    rowwise() %>%
    mutate(termgene = ifelse(is.na(termgene), "", termgene))
  # annot.out$terms.annot$topic <- sapply(annot.out$terms.annot$topic, function(x) paste("topic_", x, sep = ""))
  annot.out$terms.annot$gene <- sapply(annot.out$terms.annot$termgene, function(x){
    jsplit <- strsplit(x, ";")[[1]]
    if (length(jsplit) > 1){
      return(jsplit[[2]])
    } else {
      return("")
    }
  })
  
  annot.out$terms.annot <- annot.out$terms.annot %>% 
    group_by(topic) %>%
    mutate(rnk.rev = rank(weight))
  
  topics.mat.tmp <- left_join(dat.umap, data.frame(cell = rownames(tm.result.glm$topics), tm.result.glm$topics, stringsAsFactors = FALSE))
  
  outpdf <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/zebrafish_scchic/celltyping/from_GLMPCA/celltyping_plots_GLMPCA_", jmark, ".pdf")
  pdf(outpdf, width = 1020/72, height = 815/72, useDingbats = FALSE)
  for (i in seq(ncol(tm.result.glm$topics))){
    jtop <- paste0("dim", i)
    print(jtop)
    # plot UMAP
    
    m.umap <- PlotXYWithColor(topics.mat.tmp, xvar = "umap1", yvar = "umap2", cname = jtop, cont.color = TRUE, jtitle = jtop) + 
      scale_color_viridis_c()
    m.umap.rev <- PlotXYWithColor(topics.mat.tmp, xvar = "umap1", yvar = "umap2", cname = jtop, cont.color = TRUE, jtitle = jtop) + 
      scale_color_viridis_c(direction = -1)
    
    jsub.terms <- subset(annot.out$terms.annot, topic == jtop) %>%
      dplyr::filter(rnk <= keepn)
    
    jsub.terms.rev <- subset(annot.out$terms.annot, topic == jtop) %>%
      dplyr::filter(rnk.rev <= keepn)
    
    m.top <- jsub.terms %>%
      ungroup() %>%
      mutate(term = forcats::fct_reorder(term, dplyr::desc(weight))) %>%
      ggplot(aes(x = term, y = weight, label = gene)) +
      geom_point(size = 0.25) +
      theme_bw(8) +
      geom_text_repel(size = 2.5, segment.size = 0.1, segment.alpha = 0.25) +
      theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 2)) +
      xlab("") + ylab("Log10 Bin Weight") +
      ggtitle(paste("Top peak weights for:", jtop))
    
    m.top.rev <- jsub.terms.rev %>%
      ungroup() %>%
      mutate(term = forcats::fct_reorder(term, weight)) %>%
      ggplot(aes(x = term, y = weight, label = gene)) +
      geom_point(size = 0.25) +
      theme_bw(8) +
      geom_text_repel(size = 2.5, segment.size = 0.1, segment.alpha = 0.25) +
      theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 2)) +
      xlab("") + ylab("Log10 Bin Weight") +
      ggtitle(paste("Top peak weights for:", jtop))
    
    # check gene expression across genes
    gfilt <- unique(jsub.terms$gene)
    gfilt.rev <- unique(jsub.terms.rev$gene)
    # m.ctype <- subset(dat.bulk, gene %in% gfilt) %>%
    
    dat.bulk.sub <- subset(dat.bulk, gene %in% gfilt)
    dat.bulk.sub.rev <- subset(dat.bulk, gene %in% gfilt.rev)
    ngenes.kept <- length(unique(dat.bulk.sub$gene))
    ngenes.kept.rev <- length(unique(dat.bulk.sub.rev$gene))
    m.ctype <- dat.bulk.sub %>%
      ggplot(., aes(x = forcats::fct_reorder(celltype, dplyr::desc(zscore), .fun = median), y = zscore)) + geom_boxplot() + 
      geom_jitter(width = 0.2) + 
      theme_bw(24) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
      ggtitle(jtop, paste("Ngenes:", ngenes.kept))
    m.ctype.rev <- dat.bulk.sub.rev %>%
      ggplot(., aes(x = forcats::fct_reorder(celltype, dplyr::desc(zscore), .fun = median), y = zscore)) + geom_boxplot() + 
      geom_jitter(width = 0.2) + 
      theme_bw(24) + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + 
      ggtitle(jtop, paste("Ngenes:", ngenes.kept.rev))
    
    print(m.umap)
    print(m.umap.rev)
    print(m.top)
    print(m.top.rev)
    print(m.ctype)
    print(m.ctype.rev)
    
  }
  dev.off()
  
  
  
}


