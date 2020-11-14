# Jake Yeung
# Date of Creation: 2020-11-05
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/4-LDA_downstream_cleaner.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)

library(topicmodels)

library(hash)
library(igraph)
library(umap)

library(JFuncs)
library(scchicFuncs)

library(ggrepel)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

# Load LDA output ---------------------------------------------------------

# cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597", "#a4490d", "#ba0989")
cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf", "#28f9ff", "#88497e", "#bcf5c3", "#86f115", "#c3c89d", "#ff010b", "#664754", "#2af022", "#3afde0", "#b9b2a8", "#f6af7c", "#c3f582", "#3b3a9e", "#71a1ee", "#df5ba4", "#3a592e", "#010233", "#686cc2", "#9b114d", "#e6e6ba", "#b9f6c5")
jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

# jsuffix <- "from_hiddendomains"
jsuffix <- "from_TES"
# indir <- "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.clusterfilt.2020-11-04"
indir <- paste0("jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_BMround2all_MergeWithOld.", jsuffix)


outdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/LDA_downstream.BMround2all.merge_with_old.cleaner.", jsuffix)
dir.create(outdir)

outrdata <- file.path(outdir, paste0("LDA_downstream_objects.", Sys.Date(), ".again.", jsuffix, ".RData"))
assertthat::assert_that(!file.exists(outrdata))


# lda_outputs.count_mat_from_hiddendomains.H3K4me3.K-30.binarize.FALSE/
infs <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- file.path(hubprefix, indir, paste0("lda_outputs.count_mat_", jsuffix, ".", jmark, ".K-30.binarize.FALSE/ldaOut.count_mat_", jsuffix, ".", jmark, ".K-30.Robj"))
  print(inf)
  assertthat::assert_that(file.exists(inf))
  return(inf)
})

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

lda.outs <- lapply(infs, function(inf){
  print(inf)
  load(inf, v=T)  # out.lda
  tm.result.tmp <- posterior(out.lda)
  tm.result.tmp <- AddTopicToTmResult(tm.result.tmp, jsep = "")
  # add chr to terms
  colnames(tm.result.tmp$terms) <- paste("chr", colnames(tm.result.tmp$terms), sep = "")
  # remove suffix
  # colnames(tm.result.tmp$terms) <- sapply(colnames(tm.result.tmp$terms), function(x) strsplit(x, ";")[[1]][[1]])
  dat.umap.tmp <- DoUmapAndLouvain(tm.result.tmp$topics, jsettings)
  return(list(tm.result = tm.result.tmp, dat.umap = dat.umap.tmp))
})

dats.umap <- lapply(lda.outs, function(jout){
  jout$dat.umap
})

tm.result.lst <- lapply(lda.outs, function(jout){
  jout$tm.result
})

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")

# dat.vars.merge <- lapply(tm.result.lst, function(tm.result){
#   dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
#   # remove duplicate names
#   good.rows <- !duplicated(rownames(dat.impute.log))
#   dat.impute.log <- dat.impute.log[good.rows, ]
#   # make rownames unique 
#   # rownames(dat.impute.log) <- names(rownames(dat.impute.log))
#   dat.var <- CalculateVarAll(dat.impute.log, jchromos)
#   return(dat.var)
# }) %>%
#   bind_rows()

# dat.vars.merge <- lapply(tm.result.lst, function(tm.result){
#   dat.impute.log <- log2(t(tm.result$topics %*% tm.result$terms))
#   dat.var <- CalculateVarAll(dat.impute.log, jchromos)
#   return(dat.var)
# }) %>%
#   bind_rows()

# Plot outputs ------------------------------------------------------------

m.lst <- lapply(jmarks, function(jmark){
  jdat <- dats.umap[[jmark]]
  m <- ggplot(jdat, aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    scale_color_manual(values = cbPalette) + 
    ggtitle(jmark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m)
})

# label celltypes and replot 
# "/hpc/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs.merge_with_new_BM/clusterfilt.2020-11-04/cell_cluster_table.old_merged_with_new.H3K4me1.remove_bad_clusters.2020-11-04.txt"

indir.annot <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/count_mat_all_and_HSCs.merge_with_new_BM/clusterfilt.2020-11-04")
dats.annot <- lapply(jmarks, function(jmark){
  fname <- paste0("cell_cluster_table.old_merged_with_new.", jmark, ".remove_bad_clusters.2020-11-04.txt")
  inf.annot <- file.path(indir.annot, fname)
  dat.annot <- fread(inf.annot)
  return(dat.annot)
})

dats.cellannot <- lapply(dats.annot, function(jdat) jdat %>% dplyr::select(c(cell, cluster)))

dats.umap.annot <- lapply(jmarks, function(jmark){
  left_join(dats.umap[[jmark]], dats.cellannot[[jmark]]) %>%
    rowwise() %>%
    mutate(plate = ClipLast(cell, jsep = "_"),
           experi = ClipLast(cell, jsep = "-"),
           mark = jmark)
})


m.lst <- lapply(jmarks, function(jmark){
  jdat <- dats.umap.annot[[jmark]]
  m <- ggplot(jdat, aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + 
    scale_color_manual(values = cbPalette) + 
    ggtitle(jmark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m)
})

# Check batch efffects ----------------------------------------------------

m.batches.lst <- lapply(jmarks, function(jmark){
  jdat <- dats.umap.annot[[jmark]]
  m <- ggplot(jdat, aes(x = umap1, y = umap2, color = cluster)) + 
    geom_point() + 
    facet_wrap(~experi) + 
    scale_color_manual(values = cbPalette) + 
    ggtitle(jmark) + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  return(m)
})



# Do LDA downstream stuff -------------------------------------------------

dat.merge <- dats.umap.annot %>%
  bind_rows()
  # left_join(dat.vars.merge, by = "cell")

print("dat.merge")
print(head(dat.merge))

# Annotate from Giladi ----------------------------------------------------


jinf.tss <- file.path(hubprefix, "jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed")
assertthat::assert_that(file.exists(jinf.tss))

inf.annot <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Giladi_et_al_2018/giladi_pseudobulk_datsum_and_ncells.DESeq_and_QuantNorm.RData")
load(inf.annot, v=T)

dat.sum.long <- data.table::melt(dat.sum.norm.quantnorm)
colnames(dat.sum.long) <- c("gene", "celltype", "exprs")

dat.sum.long <- dat.sum.long %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE)) %>%
  filter(!is.nan(zscore))

topn <- 150

for (jmark in jmarks){
  
  outpdf <- file.path(outdir, paste0("LDA_downstream_celltypes_Giladi.", jmark, ".again2.pdf"))
  
  tm.result <- tm.result.lst[[jmark]]
  topics.mat <- tm.result$topics
  topics.sum <- OrderTopicsByEntropy(tm.result, jquantile = 0.99)
  annots.out <- AnnotateBins2(tm.result$terms, top.thres = 0.995, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db")
  annots.out$terms.annot$gene <- sapply(annots.out$terms.annot$termgene, function(x) ifelse(is.na(x), NA, strsplit(x, ";")[[1]][[2]]))
  
  dat.topics <- data.frame(cell = rownames(topics.mat), topics.mat, stringsAsFactors = FALSE)
  
  dat.umap.long.merge <- left_join(dat.merge %>% filter(mark == jmark), dat.topics, by = "cell")
  
  m.plates <- ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = as.character(plate))) + 
    geom_point() + 
    theme_bw() + 
    scale_color_manual(values = cbPalette) + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m.louvains <- ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = as.character(louvain))) + 
    geom_point() + 
    theme_bw() + 
    scale_color_manual(values = cbPalette) + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # m.stypes <- ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = as.character(stype))) + 
  #   geom_point() + 
  #   theme_bw() + 
  #   scale_color_manual(values = cbPalette) + 
  #   ggtitle(jmark) + 
  #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # m.var <- ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
  #   geom_point() + 
  #   theme_bw() + 
  #   scale_color_viridis_c(direction = -1) + 
  #   ggtitle(jmark) + 
  #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # m.l2r <- ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = log2(chromocounts / spikeincounts))) + 
  #   geom_point() + 
  #   theme_bw() + 
  #   scale_color_viridis_c(direction = 1) + 
  #   ggtitle(jmark) + 
  #   theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  
  # Plot a topikkukkc ------------------------------------------------------------
  
  jtopic <- topics.sum$topic[[1]]
  
  pdf(outpdf, useDingbats = FALSE)
  
  print(m.plates)
  print(m.plates + facet_wrap(~experi))
  print(m.louvains)
  print(m.louvains + facet_wrap(~experi))
  # print(m.stypes)
  # print(m.stypes + facet_wrap(~plate))
  # print(m.var)
  # print(m.var + facet_wrap(~experi))
  # print(m.l2r)
  # print(m.l2r + facet_wrap(~plate))
  
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

# aveoutputs
# save(dat.merge, dat.spikeins.mat, tm.result.lst, file = outrdata)
save(dat.merge, tm.result.lst, file = outrdata)






