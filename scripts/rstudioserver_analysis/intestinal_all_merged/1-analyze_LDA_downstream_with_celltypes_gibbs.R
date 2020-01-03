# Jake Yeung
# Date of Creation: 2019-12-31
# File: ~/projects/scchic/scripts/rstudioserver_analysis/intestinal_all_merged/1-analyze_LDA_downstream_with_celltypes_gibbs.R
# 

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

library(hash)
library(igraph)
library(umap)

library(ggrepel)

library(scchicFuncs)

library(topicmodels)

jstart <- Sys.time()

# Load LDA ----------------------------------------------------------------
ncores <- 20

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

keeptop <- 150
inf.tss <- "/home/jyeung/hpc/databases/gene_tss/gene_tss_winsize.50000.bed"
# jmark <- "k4me1"
jmethod <- "gibbs"
# assertthat::assert_that(jmethod == "vi" | jmethod == "gibbs")

inmain <- "/home/jyeung/hpc/intestinal_scchic/LDA_outputs/topicmodels/ldaAnalysisBins_intestines.2019-12-22"
infs <- list.files(inmain, pattern = ".Robj", recursive = TRUE, all.files = TRUE, full.names = TRUE)

outmain <- paste0("/home/jyeung/hpc/intestinal_scchic/from_rstudioiserver/pdfs/LDA_downstream_gibbs_withLouvain2")
dir.create(outmain)

cbPalette <- c("#696969", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#0070ea", "#15ef29", "#e1071e", "#870f86", "#cb255d")

# load pseudobulk
inf.bulk <- "/home/jyeung/hpc/intestinal_scchic/public_data/pseudobulk_data/pseudobulk_Haber_et_al_intestinal_celltypes_2019-12-30.VST_noQuantNorm.rds"
dat.mat <- readRDS(inf.bulk) %>%
as.data.frame()
dat.mat$gene <- rownames(dat.mat)

# quaneil normalize?
# library(preprocessCore)
dat.norm.long <- tidyr::gather(dat.mat, key = "celltype", value = "exprs", -gene)
# colnames(dat.norm.long) <- c("gene", "celltype", "VSTexprs")
dat.norm.long <- dat.norm.long %>%
group_by(gene) %>%
mutate(zscore = scale(exprs, center = TRUE, scale = TRUE))

# for (inf in infs){
parallel::mclapply(infs, function(inf){

  fbase <- gsub("\\", "", ClipLast(basename(inf), jsep = "\\."), fixed = TRUE)
  dbase <- dirname(inf)
  jbin <- strsplit(dbase, "\\.")[[1]][length(strsplit(dbase, "\\.")[[1]])]
  
  # outpdf <- paste0("/home/jyeung/hpc/intestinal_scchic/from_rstudioiserver/pdfs/LDA_downstream/intestines_", prefix, ".", jmark, ".", jmethod, ".pdf")
  outpdf <- file.path(outmain, paste0(fbase, ".binarize_", jbin, ".pdf"))
  if (file.exists(outpdf)){
    print(paste("Exists, skipping, ", outpdf))
    next
  }
  
  print(paste("Summarizing", inf))
  load(inf, v=T)
  tm.result <- posterior(out.lda)
  colnames(tm.result$topics) <- paste("topic", colnames(tm.result$topics), sep = "_")
  
  # Run UMAP  ---------------------------------------------------------------
  
  umap.out <- umap(tm.result$topics, config = jsettings)
  dat.umap.long <- data.frame(cell = rownames(umap.out$layout), umap1 = umap.out$layout[, 1], umap2 = umap.out$layout[, 2], stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(is.stem = grepl(pattern = "stemcells", x = cell),
           is.lgr5 = grepl(pattern = "Lgr5", x = cell),
           is.stem = is.stem | is.lgr5,
           experi = ClipLast(cell, jsep = "_"),
           prefix = paste(strsplit(cell, "-")[[1]][1:3], collapse = "-"))
  unique(dat.umap.long$experi)
  
  # add variance
  dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
  jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
  dat.var <- CalculateVarAll(dat.impute.log, jchromos)
  dat.umap.long <- left_join(dat.umap.long, dat.var)
  
  
  # Get bns -----------------------------------------------------------------
  
  annot.out <- AnnotateBins(terms.mat = tm.result$terms, inf.tss = inf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db")
  
  terms.sum <- annot.out$terms.filt %>%
    group_by(gene) %>%
    dplyr::filter(rnk == min(rnk))
  
  term2gene <- hash(terms.sum$term, terms.sum$gene)
  gene2term <- hash(terms.sum$gene, terms.sum$term)
  
  terms.filt <- annot.out$terms.filt %>%
    mutate(topic = paste("topic_", topic, sep = ""))
  
  
  # Find celltypes  ---------------------------------------------------------
  
  
  # order topics by entropy 
  topics.sum <- OrderTopicsByEntropy(tm.result = tm.result, jquantile = 0.95)
  
  # plot a topic 
  dat.merged.topics <- left_join(dat.umap.long, data.frame(cell = rownames(tm.result$topics), tm.result$topics, stringsAsFactors = FALSE))
  
  jtop <- topics.sum$topic[[1]]
  PlotXYWithColor(dat.merged.topics, xvar = "umap1", yvar = "umap2", cname = jtop)
  
  # get top genes associated with topic
  
  print(jtop)
  
  pdf(outpdf, useDingbats = FALSE)
  
  m.prefix <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = prefix)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  m.prefix.facet <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = prefix)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~prefix)
  m.experi <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  m.experi.facet <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = experi)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~prefix)
  m.var <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c()
  m.var.facet <- ggplot(dat.umap.long, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    facet_wrap(~prefix) + 
    scale_color_viridis_c()
  print(m.prefix)
  print(m.prefix.facet)
  print(m.experi)
  print(m.experi.facet)
  print(m.var)
  print(m.var.facet)
  
  for (jtop in topics.sum$topic){
    m.umap <- PlotXYWithColor(dat.merged.topics, xvar = "umap1", yvar = "umap2", cname = jtop)
    
    top.genes <- subset(terms.filt, topic == jtop & rnk <= keeptop)$gene
    assertthat::assert_that(length(top.genes) > 0)
    
    jsub <- subset(dat.norm.long, gene %in% top.genes)
    jsub.sorted.summarised <- jsub %>% group_by(celltype) %>% 
      summarise(zscore = median(zscore)) %>% 
      arrange(desc(zscore)) %>% dplyr::select(celltype)
    jlevels <- as.character(jsub.sorted.summarised$celltype)
    jsub$celltype <- factor(jsub$celltype, levels = jlevels)
    m.exprs <- ggplot(jsub,
                      aes(x = celltype , y = zscore)) +
      geom_boxplot(outlier.shape = NA) +
      # geom_violin() +
      geom_jitter(width = 0.1, size = 0.5) +
      # geom_line() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
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
print(Sys.time() - jstart)
}, 
mc.cores = ncores)

