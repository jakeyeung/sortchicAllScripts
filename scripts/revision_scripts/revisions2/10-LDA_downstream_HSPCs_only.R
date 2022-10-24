# Jake Yeung
# Date of Creation: 2022-07-29
# File: ~/projects/scchic/scripts/revision_scripts/revisions2/10-LDA_downstream_HSPCs_only.R
# 

rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(topicmodels)
library(scchicFuncs)

library(ggrepel)

library(ggrastr)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)

library(hash)
library(igraph)
library(umap)

jsettings <- umap.defaults
jsettings[["n_neighbors"]] <- 30
jsettings[["min_dist"]] <- 0.1
jsettings[["random_state"]] <- 123



AnnotateBins2.R4 <- function (terms.mat, top.thres = 0.995, inf.tss = "/Users/yeung/data/scchic/tables/gene_tss_winsize.50000.bed", 
                              txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", 
                              chromos.keep = c(paste("chr", seq(19), sep = ""), "chrX", 
                                               "chrY"), skip.split = FALSE) 
{
  assertthat::assert_that(file.exists(inf.tss))
  # assertthat::assert_that(class(terms.mat) == "matrix")
  regions <- data.frame(seqnames = sapply(colnames(terms.mat), 
                                          GetChromo), start = sapply(colnames(terms.mat), GetStart), 
                        end = sapply(colnames(terms.mat), GetEnd), stringsAsFactors = FALSE)
  rownames(regions) <- colnames(terms.mat)
  print("Chromos to keep")
  print(chromos.keep)
  regions <- subset(regions, seqnames %in% chromos.keep)
  regions.range <- makeGRangesFromDataFrame(as.data.frame(regions))
  regions.annotated <- as.data.frame(annotatePeak(regions.range, 
                                                  TxDb = txdb, annoDb = annodb))
  regions.annotated$region_coord <- names(regions.range)
  topic.regions <- lapply(seq(nrow(terms.mat)), function(clst) {
    return(SelectTopRegions(terms.mat[clst, ], colnames(terms.mat), 
                            method = "thres", method.val = top.thres))
  })
  print(paste("Using TSS definitions from:", inf.tss))
  terms.long <- data.frame(term = colnames(terms.mat), as.data.frame(t(terms.mat)), 
                           stringsAsFactors = FALSE) %>% gather(key = "topic", value = "weight", 
                                                                -term) %>% mutate(topic = gsub("X", "", topic)) %>% group_by(topic) %>% 
    arrange(desc(weight)) %>% mutate(rnk = seq(length(weight))) %>% 
    rowwise()
  terms.filt.top <- terms.long %>% rowwise()
  tss.dat <- fread(inf.tss, col.names = c("seqnames", "start", 
                                          "end", "tssname"))
  tss.dat$gene <- sapply(tss.dat$tssname, function(x) strsplit(x, 
                                                               ";")[[1]][[2]])
  annots.biomart <- regions.annotated %>% mutate(midpt = start + 
                                                   (end - start)/2) %>% filter(region_coord %in% terms.filt.top$term)
  annots.gr <- makeGRangesFromDataFrame(annots.biomart %>% 
                                          dplyr::select(seqnames, start, end, SYMBOL, region_coord, 
                                                        ENSEMBL), keep.extra.columns = TRUE)
  annots.tss.gr <- makeGRangesFromDataFrame(tss.dat, keep.extra.columns = TRUE)
  out <- findOverlaps(annots.tss.gr, annots.gr, type = "within")
  out2 <- findOverlaps(annots.gr, annots.tss.gr, type = "any")
  out2.df = data.frame(annots.gr[queryHits(out2), ], annots.tss.gr[subjectHits(out2), 
  ]) %>% mutate(midpt = start + round(width/2), midpt.1 = start.1 + 
                  round(width.1/2), dist.to.tss = midpt.1 - midpt)
  out2.df.closest <- out2.df %>% group_by(region_coord) %>% 
    filter(abs(dist.to.tss) == min(abs(dist.to.tss)))
  terms.new <- paste(out2.df.closest$region_coord, out2.df.closest$gene, 
                     sep = ";")
  terms.hash <- hash::hash(out2.df.closest$region_coord, terms.new)
  terms.annot <- terms.filt.top %>% mutate(termgene = ifelse(!is.null(terms.hash[[term]]), 
                                                             terms.hash[[term]], NA))
  terms.filt <- terms.filt.top %>% mutate(termgene = ifelse(!is.null(terms.hash[[term]]), 
                                                            terms.hash[[term]], NA)) %>% filter(!is.na(termgene))
  if (!skip.split) {
    terms.filt <- terms.filt %>% mutate(gene = sapply(termgene, 
                                                      function(x) strsplit(x, ";")[[1]][[2]])) %>% group_by(gene)
  }
  return(list(topic.regions = topic.regions, regions.annotated = regions.annotated, 
              terms.annot = terms.annot, out2.df.closest = out2.df.closest, 
              terms.filt = terms.filt))
}



# jmarks <- c("k4me1", "k4me3", "k27me3", "k9me3"); names(jmarks) <- jmarks
# jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarksold) <- jmarks
jmarks <- c("k4me1", "k4me3", "k27me3"); names(jmarks) <- jmarks
jmarksold <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarksold) <- jmarks

jinf.tss <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databases/gene_tss/gene_tss_winsize.50000.bed")
assertthat::assert_that(file.exists(jinf.tss))

# Load new colors ---------------------------------------------------------

inf.colors.fixed <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/ctypes_on_umap_batch_corrected_colors_fixed/dat_colors_DC_monocyte_fixed.2022-05-17.txt"
dat.colors.fixed <- fread(inf.colors.fixed)


# Load meta ----------------------------------------------------------------

dat.meta.lst <- lapply(jmarks, function(jmark){
  inf.meta <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/from_jupyterhub/primetime_plots/umaps_pcas_with_batch_corrections/umap_metadata_primetime.", jmark, ".2022-04-21.txt")
  dat.meta <- fread(inf.meta) %>%
    left_join(., dat.colors.fixed) %>%
    rowwise() %>%
    mutate(colcode = colcodenew)
  # replace colcode with colcodenew
})

dat.meta.colors <- subset(dat.meta.lst$k4me1, select = c(ctype.from.LL, colcode))
ctype2col <- hash::hash(dat.meta.colors$ctype.from.LL, dat.meta.colors$colcode)




# Load LDA  ---------------------------------------------------------------

lda.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_progs_only/lda_outputs.countmat_TSS_new_progenitors_only.", jmark, ".2022-07-27/ldaOut.countmat_TSS_new_progenitors_only.", jmark, ".2022-07-27.Robj")
  load(inf, v=T)
  return(out.lda)
})

count.mat.norm.lst <- lapply(jmarks, function(jmark){
  print(jmark)
  inf <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/LDA_outputs_progs_only/lda_outputs.countmat_TSS_new_progenitors_only.", jmark, ".2022-07-27/ldaOut.countmat_TSS_new_progenitors_only.", jmark, ".2022-07-27.Robj")
  load(inf, v=T)
  count.mat.norm <- sweep(count.mat, MARGIN = 2, STATS = colSums(count.mat), FUN = "/")
  return(count.mat.norm)
})

tm.lst <- lapply(lda.lst, function(lda){
  posterior(lda) %>%
    AddTopicToTmResult()
})


dat.meta.lst <- lapply(jmarks, function(jmark){
  dat.umap <- DoUmapAndLouvain(tm.lst[[jmark]]$topics, jsettings = jsettings) %>%
    rowwise() %>%
    mutate(platename = ClipLast(cell, jsep = "_")) %>%
    left_join(., dat.meta.lst[[jmark]] %>% dplyr::select(c(-umap1, -umap2, -louvain)))
})

m.lst <- lapply(jmarks, function(jmark){
  ggplot(dat.meta.lst[[jmark]], aes(x = umap1, y = umap2, color = colcode)) + 
    geom_point() + 
    theme_bw() + 
    ggtitle(jmark) + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode, guide = "legend") +
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
})

# check k27me3 effects
m.lst$k27me3 + facet_wrap(~platename)

# Plot topic loadings -----------------------------------------------------

outdir <- "/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/scChiC/new_experiments/downstream_LDA_progs_only"

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")


inf.annot <- paste0("/nfs/scistore12/hpcgrp/jyeung/data_from_Hubrecht/hpc_hub_oudenaarden/databases/celltyping_reference_data/giladi_pseudobulk_datsum_and_ncells.DESeq_and_QuantNorm.WithC1qb.2021-11-02.RData")
load(inf.annot, v=T)

dat.sum.long <- data.table::melt(dat.sum.norm.quantnorm)
colnames(dat.sum.long) <- c("gene", "celltype", "exprs")

dat.sum.long <- dat.sum.long %>%
  group_by(gene) %>%
  mutate(zscore = scale(exprs, center = TRUE, scale = TRUE)) %>%
  filter(!is.nan(zscore))

topn <- 150
for (jsuffix in jmarks){
  print(jsuffix)
  outpdf <- file.path(outdir, paste0("LDA_downstream_celltypes_Giladi.", jsuffix, ".pdf"))
  
  tm.result <- tm.lst[[jsuffix]]
  # tm.result <- tm.lst[[jsuffix]]
  # tm.result <- AddTopicToTmResult(tm.result)
  topics.mat <- tm.result$topics
  topics.sum <- OrderTopicsByEntropy(tm.result, jquantile = 0.99)
  terms.mat.clean <- tm.result$terms
  print(dim(terms.mat.clean))
  colnames(terms.mat.clean) <- paste("chr", sapply(colnames(terms.mat.clean), function(x) strsplit(x, ";")[[1]][[1]]), sep = "")
  # remove dupes
  terms.mat.clean <- terms.mat.clean[, !duplicated(colnames(terms.mat.clean))]
  print(dim(terms.mat.clean))
  annots.out <- AnnotateBins2.R4(terms.mat.clean, top.thres = 0.995, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db")
  annots.out$terms.annot$gene <- sapply(annots.out$terms.annot$termgene, function(x) ifelse(is.na(x), NA, strsplit(x, ";")[[1]][[2]]))
  
  dat.topics <- data.frame(cell = rownames(topics.mat), topics.mat, stringsAsFactors = FALSE)
  
  dat.umap.long.merge <- left_join(dat.meta.lst[[jsuffix]], dat.topics, by = "cell") 
  
  m.plates <- ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = colcode)) + 
    geom_point() + 
    theme_bw() + 
    scale_color_identity( labels = dat.meta.colors$ctype.from.LL, breaks = dat.meta.colors$colcode, guide = "legend") +
    ggtitle(jsuffix) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  m.louvains <- ggplot(dat.umap.long.merge, aes(x = umap1, y = umap2, color = as.character(louvain))) + 
    geom_point() + 
    theme_bw() + 
    scale_color_manual(values = cbPalette) + 
    ggtitle(jsuffix) + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  # Plot a topikkukkc ------------------------------------------------------------
  
  jtopic <- topics.sum$topic[[1]]
  
  pdf(outpdf, useDingbats = FALSE)
  
  print(m.lst[[jsuffix]])
  # dprint(m.var.lst[[jsuffix]])
  
  print(m.plates)
  print(m.plates + facet_wrap(~ctype.from.LL))
  print(m.louvains)
  print(m.louvains + facet_wrap(~ctype.from.LL))
  
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
      geom_text_repel(size = 2, segment.size = 0.1, segment.alpha = 0.25, max.overlaps = 50) +
      theme(aspect.ratio=0.2, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 3.5)) +
      xlab("") + ylab("Log10 Bin Weight") +
      ggtitle(paste("Top peak weights for:", jtopic))
    print(m.top)
  }
  dev.off()
}


# Write umaps  ------------------------------------------------------------

for (jmark in jmarks){
  print(jmark)
  outftmp <- file.path(outdir, paste0("dat_umap_progs_only.", jmark, ".", Sys.Date(), ".txt"))
  fwrite(dat.meta.lst[[jmark]], file = outftmp, sep = "\t")
}


# Look at marker genes ----------------------------------------------------


marker.genes.c1 <- c("Tgm2", "Mmrn1", "Trim47", "Pdzk1ip1", "Mllt3", "Mecom", "Esam", "Cish", "Hes1", "Mycn")
marker.genes.c2 <- c("Dntt", "Gm5111", "Flt3", "9030619P08Rik", "Il1r1", "Vldlr")
marker.genes.c3 <- c("Uhrf1", "Rmr2", "Lig1", "Tipin", "Rad54l", "Fabp5", "Nrm", "Gatm", "Hn1l", "Fam111a")
marker.genes.mk <- c("Pf4", "Gata1", "Slc14a1", "F2r", "Itga2b", "Zfp385a", "Zfpm1", "Plek", "Cd9", "Zeb2")

marker.genes <- list(marker.genes.c1, marker.genes.c2, marker.genes.c3, marker.genes.mk)
names(marker.genes) <- c("c1", "c2", "c3", "mk")


# Plot marker  genes  -----------------------------------------------------


library(DescTools)
jnames <- names(marker.genes)

# jname <- "c2"
for (markref in jmarks){
  outpdf <- file.path(outdir, paste0("marker_genes_C1_C2_C3.", markref, ".", Sys.Date(), ".pdf"))
  pdf(outpdf, useDingbats = FALSE)
  for (jname in jnames){
    genes.grep <- paste0(marker.genes[[jname]], collapse = "|")
    count.mat <- count.mat.norm.lst[[markref]]
    jsignal <- colMeans(count.mat[grepl(genes.grep, rownames(count.mat)), ])
    # winsorize
    jsignal.win <- DescTools::Winsorize(jsignal, probs = c(0.05, 0.95))
    dat.signal <- data.frame(cell = names(jsignal), frac.counts = jsignal, frac.counts.win = jsignal.win, stringsAsFactors = FALSE) %>%
      left_join(., dat.meta.lst[[markref]])
    
    m <- ggplot(dat.signal, aes(x = umap1, y = umap2, color = log2(frac.counts.win * 1000 + 1))) + 
      geom_point() + 
      ggtitle(markref, jname) + 
      scale_color_viridis_c() + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
    print(m)
  }
  dev.off()
}

