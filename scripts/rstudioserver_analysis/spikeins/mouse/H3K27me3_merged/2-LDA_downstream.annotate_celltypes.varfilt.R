# Jake Yeung
# Date of Creation: 2020-10-04
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/mouse/H3K27me3_merged/2-LDA_downstream.annotate_celltypes.R
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
library(scchicFuncs)

library(topicmodels)

library(ggrepel)

# Functions ---------------------------------------------------------------


AnnotateSortFromLayout <- function(plate, rowcoord, colcoord){
  assertthat::assert_that(is.numeric(plate))
  assertthat::assert_that(is.numeric(rowcoord))
  assertthat::assert_that(is.numeric(colcoord))
  if (rowcoord >= 1 & rowcoord <= 8 & colcoord == 1){
    print("Empty cell")
    warning("Empty cell")
    ctype <- "Empty"
    return(ctype)
  }
  if (plate >= 1 & plate <= 7){
    if (colcoord >= 1 & colcoord <= 11){
      ctype <- "Unenriched"
    } else if (colcoord >= 12 & colcoord <= 18){
      ctype <- "LinNeg"
    } else if (colcoord >= 19 & colcoord <= 24){
      ctype <- "LSK"
    }
    
  } else if (plate >= 8 & plate <= 13){
    if (colcoord >= 1 & colcoord <= 12){
      ctype <- "Unenriched"
    } else if (colcoord >= 13 & colcoord <= 24){
      ctype <- "LinNeg"
    }
  } else {
    print(paste("Unknown plate:", plate))
    ctype <- "Unknown"
  }
  return(ctype)
}



# Constnats ---------------------------------------------------------------

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)



cbPalette <- c("#696969", "#32CD32", "#56B4E9", "#FFB6C1", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400", "#FFB6C1", "#32CD32", "#0b1b7f", "#ff9f7d", "#eb9d01", "#7fbedf")

# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks
jmarks <- c("H3K27me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"

# inmain <- file.path(hubprefix, "jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_VAN5046")
# inmain <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_VAN5046_varfilt"
# inmain <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_VAN5230_BM"
# 
# inmain.suffix <- "VAN5046_VAN5230_BM"
# inmain.suffix <- "VAN5230_BM"
inmain.suffix <- "VAN5046_VAN5230_BM_varfilt"
inmain <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/raw_demultiplexed/LDA_outputs_all_spikeins/ldaAnalysisBins_mouse_spikein_", inmain.suffix)
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
topn <- 150

# jvarfilt <- 1
# jsuffix <- "_counts_filt.2020-09-12.K-30"
# jsuffix <- file.path(paste0("lda_outputs.count_mat_", jmark, "_l2r_filt.2020-09-13.minl2r_-1.varfilt_", jvarfilt, ".K-30.binarize.FALSE"), paste0("ldaOut.count_mat_", jmark, "_l2r_filt.2020-09-13.minl2r_-1.varfilt_", jvarfilt, ".K-30.Robj"))
# jsuffix <- paste0("_l2r_filt.2020-09-13.minl2r_-1.varfilt_", jvarfilt, ".K-30")
jsuffix <- paste0("_l2r_filt.2020-10-03.minl2r_0.varfilt_1.5.K-30")

outdir <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/spikeins_mouse.round2/", inmain.suffix, "/suffix", jsuffix)
dir.create(outdir, recursive = TRUE)
outpdf.check <- file.path(outdir, "plate_checks.pdf")

jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

jinf.tss <- file.path(hubprefix, "jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed")
assertthat::assert_that(file.exists(jinf.tss))

inf.annot <- file.path(hubprefix, "jyeung/data/scChiC/public_data/Giladi_et_al_2018/giladi_pseudobulk_datsum_and_ncells.DESeq_and_QuantNorm.RData")
load(inf.annot, v=T)



# Do stuff ----------------------------------------------------------------


ldas.lst <- lapply(jmarks, function(jmark){
  dname <- paste0("lda_outputs.count_mat_", jmark, jsuffix, ".binarize.FALSE")
  indir <- file.path(inmain, dname)
  fname <- paste0("ldaOut.count_mat_", jmark, jsuffix, ".Robj")
  inf.lda <- file.path(indir, fname)
  assertthat::assert_that(file.exists(inf.lda))
  print(inf.lda)
  load(inf.lda, v=T)
  tm.result <- posterior(out.lda)
  tm.result <- AddTopicToTmResult(tm.result)
  dat.umap <- DoUmapAndLouvain(tm.result$topics, jsettings)
  dat.umap$mark <- jmark
  return(list(dat.umap = dat.umap, tm.result = tm.result))
})

dat.umaps.lst <- lapply(ldas.lst, function(jlda){
  jlda$dat.umap
})

tm.result.lst <- lapply(ldas.lst, function(jlda){
  jlda$tm.result
})


# tm.result.lst <- lapply(jmarks, function(jmark){
#   dname <- paste0("lda_outputs.count_mat_", jmark, jsuffix, ".binarize.FALSE")
#   indir <- file.path(inmain, dname)
#   fname <- paste0("ldaOut.count_mat_", jmark, "_counts_filt.2020-09-12.K-30.Robj")
#   inf.lda <- file.path(indir, fname)
#   assertthat::assert_that(file.exists(inf.lda))
#   print(inf.lda)
#   load(inf.lda, v=T)
#   tm.result <- posterior(out.lda)
#   tm.result <- AddTopicToTmResult(tm.result)
#   return(tm.result)
# })


for (jmark in jmarks){
  m <- ggplot(dat.umaps.lst[[jmark]], aes(x = umap1, y = umap2, color = louvain)) + 
    geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette) + 
    ggtitle(jmark)
  print(m)
}


dat.umaps.long <- bind_rows(dat.umaps.lst) %>%
  rowwise() %>%
  mutate(experi = ClipLast(cell, jsep = "_"), 
         plate = as.numeric(strsplit(experi, "-")[[1]][[6]]),
         rowcoord = AddPlateCoordinates(cell)$rowcoord,
         colcoord = AddPlateCoordinates(cell)$colcoord,
         stype = AnnotateSortFromLayout(plate, rowcoord, colcoord))

dat.meta <- subset(dat.umaps.long, select = c(cell, experi, plate, rowcoord, colcoord, stype))

for (jmark in jmarks){
  m <- ggplot(dat.umaps.long %>% filter(mark == jmark), aes(x = umap1, y = umap2, color = stype)) + 
    geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_manual(values = cbPalette) + 
    ggtitle(jmark)
  print(m)
}


# Check annotations -------------------------------------------------------

jexperis <- unique(dat.umaps.long$experi)

pdf(outpdf.check, useDingbats = FALSE)
for (jexperi in jexperis){
  m <- ggplot(subset(dat.umaps.long, experi == jexperi), aes(y = rowcoord, x = colcoord, color = stype)) + 
    geom_point(alpha = 0.5) + 
    theme_bw() + 
    scale_color_manual(values = cbPalette) + 
    theme(aspect.ratio= 16 / 24, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_x_continuous(breaks = seq(24)) + 
    # scale_y_continuous(breaks = seq(24)) + 
    scale_y_reverse() + 
    ggtitle(jexperi)
  print(m)
}
dev.off()



# Variance stuff ----------------------------------------------------------

jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
dat.var.lst <- lapply(jmarks, function(jmark){
  tm.result <- tm.result.lst[[jmark]]
  dat.impute <- t(log2(tm.result$topics %*% tm.result$terms))
  dat.var <- CalculateVarAll(dat.impute, jchromos)
}) 

dat.var <- dat.var.lst %>%
  bind_rows()

dat.merge <- left_join(dat.umaps.long, dat.var)

for (jmark in jmarks){
  m <- ggplot(dat.merge %>% filter(mark == jmark), aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + 
    theme_bw() + theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    scale_color_viridis_c(direction = -1) + 
    ggtitle(jmark)
  print(m)
}


# Annotate from Giladi ----------------------------------------------------



dat.sum.long <- data.table::melt(dat.sum.norm.quantnorm)
colnames(dat.sum.long) <- c("gene", "celltype", "exprs")

dat.sum.long <- dat.sum.long %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE)) %>%
  filter(!is.nan(zscore))


for (jmark in jmarks){
  
  outpdf <- file.path(outdir, paste0("LDA_downstream_celltypes_Giladi.", jmark, ".pdf"))
  
  tm.result <- tm.result.lst[[jmark]]
  topics.mat <- tm.result$topics
  topics.sum <- OrderTopicsByEntropy(tm.result, jquantile = 0.99)
  annots.out <- AnnotateBins2(tm.result$terms, top.thres = 0.995, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db")
  annots.out$terms.aninot$gene <- sapply(annots.out$terms.annot$termgene, function(x) ifelse(is.na(x), NA, strsplit(x, ";")[[1]][[2]]))
  # annots.out$terms.filt$topic <- paste("topic", annots.out$terms.filt$topic, sep = "")
  # annots.out$terms.annot
  
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
  
  m.stypes <- ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = as.character(stype))) + 
    geom_point() + 
    theme_bw() + 
    scale_color_manual(values = cbPalette) + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m.var <- ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = cell.var.within.sum.norm)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_viridis_c(direction = -1) + 
    ggtitle(jmark) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # Plot a topikkukkc ------------------------------------------------------------
  
  jtopic <- topics.sum$topic[[1]]
  
  pdf(outpdf, useDingbats = FALSE)
  
  print(m.plates)
  print(m.plates + facet_wrap(~plate))
  print(m.louvains)
  print(m.louvains + facet_wrap(~plate))
  print(m.stypes)
  print(m.stypes + facet_wrap(~plate))
  print(m.var)
  print(m.var + facet_wrap(~plate))
  
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




