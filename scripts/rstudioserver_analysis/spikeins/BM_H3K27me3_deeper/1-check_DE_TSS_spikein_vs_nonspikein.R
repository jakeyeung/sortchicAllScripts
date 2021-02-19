# Jake Yeung
# Date of Creation: 2021-02-02
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_H3K27me3_deeper/1-check_DE_TSS_spikein_vs_nonspikein.R
# 

rm(list=ls())


library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(Matrix)
library(scchicFuncs)
library(JFuncs)

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomicRanges)

# jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3", "H3K9me3"); names(jmarks) <- jmarks
jmarks <- c("H3K27me3"); names(jmarks) <- jmarks

hubprefix <- "/home/jyeung/hub_oudenaarden"
jmark <- "H3K27me3"
# jmark <- "H3K4me1"
outdir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/compare_DE_analysis_spikein_vs_total"

make.plots <- FALSE

# for (jmark in jmarks){
  outpdf <- file.path(outdir, paste0("compare_DE_spikein_vs_total_normalization.", jmark, ".", Sys.Date(), ".pdf"))
  
  if (make.plots){
    pdf(outpdf, useDingbats = FALSE)
  }
  
  # Load TSS  ---------------------------------------------------------------
  # 
  # inf <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3_H3K27me3_rep2_rep3reseq/poisson_fit_TSS_10000.H3K27me3.2020-12-10.newannot2.rep2_rep3seq.RData"
  # load(inf, v=T)
  # 
  # jfits <- SummarizeParamsPvalues(jfits.lst = jfits.lst, jmark = jmark, paramname = "Cluster")
  
  # print('poisson_fit_TSS_10000.H3K4me3.2021-02-02.total.newannot2.RData')
  if (jmark == "H3K27me3"){
    inf.tss.total <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3_H3K27me3_rep2_rep3reseq/poisson_fit_TSS_10000.H3K27me3.2021-02-02.newannot2.rep2_rep3seq.with_se.RData")
    inf.tss.spikeins <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3_H3K27me3_rep2_rep3reseq/poisson_fit_TSS_10000.H3K27me3.2021-02-02.newannot2.rep2_rep3seq.with_se.spikeins.RData"
  } else {
    indir.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.spikeins.no_H3K27me3_TSS_spikein_or_total"
    inf.tss.total <- file.path(indir.tss, paste0("poisson_fit_TSS_10000.", jmark, ".2021-02-02.total.newannot2.RData"))
    inf.tss.spikeins <- file.path(indir.tss, paste0("poisson_fit_TSS_10000.", jmark, ".2021-02-02.spikeins.newannot2.RData"))
  }
  assertthat::assert_that(file.exists(inf.tss.total))
  assertthat::assert_that(file.exists(inf.tss.spikeins))
  
  
  load(inf.tss.total, v=T)
  print(head(jfits.lst[[1]]))
  jfits.tss.lst.total <- jfits.lst
  params.dat.tss.total.withpval <- SummarizeParamsPvalues(jfits.tss.lst.total, jmark = jmark, paramname = "Cluster")
  
  load(inf.tss.spikeins, v=T)
  print(head(jfits.lst[[1]]))
  jfits.tss.lst.spikeins <- jfits.lst
  params.dat.tss.spikeins.withpval <- SummarizeParamsPvalues(jfits.tss.lst.spikeins, jmark = jmark, paramname = "Cluster")
  
  
  # Load DEs  ---------------------------------------------------------------
  
  
  inf.de.spikeins <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.allmarks_spikeins/poisson_fit_bins.", jmark, ".2020-12-25.newannot2.witherrors.MoreBins.newestmeta.RData")
  assertthat::assert_that(file.exists(inf.de.spikeins))
  
  inf.de.total <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3.k4me1k9me3/poisson_fit_bins.", jmark, ".2020-12-19.newannot2.witherrors.MoreBins.newestmeta.RData")
  assertthat::assert_that(file.exists(inf.de.total))
  
  load(inf.de.spikeins, v=T)
  print(head(jfits.lst[[1]]))
  jfits.lst.spikeins <- jfits.lst
  params.dat.bins.spikeins.withpval <- SummarizeParamsPvalues(jfits.lst.spikeins, jmark = jmark, paramname = "Cluster")
  
  
  load(inf.de.total, v=T)
  print(head(jfits.lst[[1]]))
  print(head(ncuts.for.fit.mark))
  jfits.lst.total <- jfits.lst
  params.dat.bins.total.withpval <- SummarizeParamsPvalues(jfits.lst.total, jmark = jmark, paramname = "Cluster")
  
  
  # Load TSS diff exprs -----------------------------------------------------
  
  inf.mat <- file.path(hubprefix, "jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/heatmap_pdfs_and_ordered_matrices/heatmap_ordered_with_labels.H3K4me1.2021-01-08.rearranged.RData")
  load(inf.mat, v=T)
  
  rnames <- rownames(mat.adj.tmp)
  jgenes <- unique(sapply(rnames, function(x) strsplit(x, "\\.")[[1]][[4]], USE.NAMES = FALSE))
  
  coords <- sapply(rnames, function(x) strsplit(x, ";")[[1]][[1]])
  dat.bed <- GetBedFromCoords(coords, add.chr = TRUE, strip.chr = FALSE) %>%
    rowwise() %>%
    mutate(Name = paste("chr", Name, sep = ""))
  dat.bed$FullName <- rnames
  dat.bed$gene <- sapply(rnames, function(x) strsplit(x, "\\.")[[1]][[4]])
  
  jbin <- dat.bed$Name[[1]]
  
  subset(params.dat.tss.spikeins.withpval, bin == jbin)
  
  
  
  # Define ctypes -----------------------------------------------------------
  
  # filter for celltype specific, then merge the two
  params.merge.tss <- left_join(params.dat.tss.total.withpval, params.dat.tss.spikeins.withpval, by = c("bin", "param")) 
  if (jmark == "H3K27me3"){
    params.merge.tss.filt <- subset(params.merge.tss, bin %in% dat.bed$Name)
  } else {
    params.merge.tss.filt <- subset(params.merge.tss, bin %in% dat.bed$FullName)
  }
  
  params.merge.bins <- left_join(params.dat.bins.total.withpval, params.dat.bins.spikeins.withpval, by = c("bin", "param")) 
  
  # annotate bins then filter for celltype specific genes
  jinf.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed"
  jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")
  jcoords.vec <- unique(params.merge.bins$bin)
  bins.annot <- AnnotateCoordsFromList.GeneWise(coords.vec = jcoords.vec, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)
  bins.annot.filt <- subset(bins.annot$out2.df,  gene %in% dat.bed$gene)
  bins.keep <- bins.annot.filt$region_coord
  
  params.merge.bins.filt <- subset(params.merge.bins, bin %in% bins.keep)
  
  
  ctypes <- unique(params.merge.tss$param); names(ctypes) <- ctypes
  
  
  
  # Compare diff exprs for H3K27me3 spikein vs nonspimein -------------------
  
  
  print(params.merge.tss.filt)
  
  m <- ggplot(params.merge.tss.filt, aes(x = estimate.x, y = estimate.y)) + 
    geom_point() + 
    xlab("logFC norm by TSS regions") + 
    ylab("logFC norm by spikeins") + 
    ggtitle("TSS") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    facet_wrap(~param) + 
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  # look at Eryths
  
  print(head(params.merge.tss.filt))
  
  for (jparam in ctypes){
    m <- ggplot(params.merge.tss.filt %>% filter(param == jparam), aes(x = estimate.x, y = estimate.y)) + 
      geom_point(alpha = 0.25) + 
      ggtitle("TSS", jparam) + 
      xlab("logFC norm by TSS regions") + 
      ylab("logFC norm by spikeins") + 
      geom_errorbarh(mapping = aes(xmin = estimate.x - se.x, xmax = estimate.x + se.x), alpha = 0.25) + 
      geom_errorbar(mapping = aes(ymin = estimate.y - se.y, ymax = estimate.y + se.y), alpha = 0.25) + 
      geom_vline(xintercept = 0, linetype = "dotted") + 
      geom_hline(yintercept = 0, linetype = "dotted") + 
      facet_wrap(~param) + 
      coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
  }
  
  
  
  # Get genes whose signs have flipped  -------------------------------------
  
  # number of flipped genes
  
  
  jsub <- subset(params.merge.tss.filt, (estimate.x - se.x > 0 & estimate.y + se.y < 0) | (estimate.x + se.x < 0 & estimate.y - se.y > 0))
  
  flipped.genes.lst <- lapply(ctypes, function(ctype){
    print(ctype)
    jsubsub <- subset(jsub %>% filter(param == ctype), (estimate.x - se.x > 0 & estimate.y + se.y < 0) | (estimate.x + se.x < 0 & estimate.y - se.y > 0))
    flipped.genes.tss <- jsubsub$bin
    return(flipped.genes.tss)
  })
  
  flipped.genes.renamed.lst <- lapply(flipped.genes.lst, function(jcoord){
    subset(dat.bed, Name %in% jcoord)$gene
  })
  
  jgenes.eryth.flipped <- sort(unique(subset(dat.bed, Name %in% flipped.genes.lst$ClusterEryths.Estimate)$gene))
  
  inf.public <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/robjs_for_heatmaps_and_umaps/giladi_pseudobulk_exprs_data.rds"
  dat.public <- readRDS(inf.public)
  
  m2c <- MarkerToCelltype()
  dat.public$celltype2 <- sapply(as.character(dat.public$celltype), function(x) m2c[[x]])
 
  jgenes <- jgenes.eryth.flipped
  jgenes <- c("Irf8") 
  jgenes.str <- paste(jgenes, collapse = "|")
  dat.public.sub <- subset(dat.public, grepl(jgenes.str, gene))
  
  library(forcats)
  
  ggplot(dat.public.sub %>% filter(celltype %in% c("core", "Hba.a2")), 
         aes(x = forcats::fct_reorder(.f = celltype2, .x = zscore, .fun = median, .desc = TRUE), y = exprs)) + 
    ggtitle('Pseudobulk mRNA data across blood cell types', paste(jgenes, collapse = ",")) +
    geom_boxplot() + 
    geom_point() + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  
  for (jparam in ctypes){
    m <- ggplot(jsub %>% filter(param == jparam & bin %in% flipped.genes.lst[[jparam]]), aes(x = estimate.x, y = estimate.y)) + 
      geom_point(alpha = 0.25) + 
      ggtitle("TSS", jparam) + 
      xlab("logFC norm by TSS regions") + 
      ylab("logFC norm by spikeins") + 
      geom_errorbarh(mapping = aes(xmin = estimate.x - se.x, xmax = estimate.x + se.x), alpha = 0.25) + 
      geom_errorbar(mapping = aes(ymin = estimate.y - se.y, ymax = estimate.y + se.y), alpha = 0.25) + 
      geom_vline(xintercept = 0, linetype = "dotted") + 
      geom_hline(yintercept = 0, linetype = "dotted") + 
      facet_wrap(~param) + 
      coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
      theme_bw() + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    return(m)
  }
  
  # What are those genes ?  -------------------------------------------------
  
  
  subset(dat.bed, Name %in% flipped.genes.lst$ClusterEryths.Estimate)$FullName
  
  jgene <- "Irf8"
  jbin <- subset(dat.bed, grepl(jgene, FullName))$Name[[1]]
  
  subset(params.merge.tss, bin == jbin)
  
  jtmp <- subset(params.merge.tss, bin %in% flipped.genes.lst$ClusterEryths.Estimate) %>% 
    arrange(desc(estimate.x ^ 2 + estimate.y ^ 2)) %>%
    left_join(., dat.bed, by = c("bin" = "Name"))
  print(jtmp)
  
  
  m <- ggplot(params.merge.tss, aes(x = estimate.x, fill = param)) + 
    geom_density() + 
    facet_wrap(~param) + 
    xlab("logFC norm by TSS regions") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5)) + 
    ggtitle("TSS") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  m <- ggplot(params.merge.tss, aes(x = estimate.y, fill = param)) + 
    geom_density() + 
    facet_wrap(~param) + 
    xlab("logFC norm by spikeins") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5)) + 
    ggtitle("TSS") + 
    theme_bw() + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  # Check bin analysis  -----------------------------------------------------     
  
  
  
  for (jparam in ctypes){
    m <- ggplot(params.merge.bins %>% filter(param == jparam & abs(estimate.x) < 5 & abs(estimate.y) < 5), aes(x = estimate.x, y = estimate.y)) + 
      geom_point(alpha = 0.25) + 
      # geom_errorbarh(mapping = aes(xmin = estimate.x - se.x, xmax = estimate.x + se.x), alpha = 0.25) + 
      # geom_errorbar(mapping = aes(ymin = estimate.y - se.y, ymax = estimate.y + se.y), alpha = 0.25) + 
      xlab("logFC norm by total") + 
      ylab("logFC norm by spikeins") + 
      geom_vline(xintercept = 0, linetype = "dotted") + 
      geom_hline(yintercept = 0, linetype = "dotted") + 
      facet_wrap(~param) + 
      coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
      theme_bw() + 
      ggtitle("bins") + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
  }
  
  m <- ggplot(params.merge.bins, aes(x = estimate.x, y = estimate.y)) + 
    geom_point(alpha = 0.25) + 
    # geom_errorbarh(mapping = aes(xmin = estimate.x - se.x, xmax = estimate.x + se.x), alpha = 0.25) + 
    # geom_errorbar(mapping = aes(ymin = estimate.y - se.y, ymax = estimate.y + se.y), alpha = 0.25) + 
    xlab("logFC norm by total") + 
    ylab("logFC norm by spikeins") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    facet_wrap(~param) + 
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
    theme_bw() + 
    ggtitle("bins") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  m <- ggplot(params.merge.bins, aes(x = estimate.x, fill = param)) + 
    geom_density() + 
    xlab("logFC norm by total") + 
    facet_wrap(~param) + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5)) + 
    theme_bw() + 
    ggtitle("bins") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  m <- ggplot(params.merge.bins, aes(x = estimate.y, fill = param)) + 
    geom_density() + 
    xlab("logFC norm by spikeins") + 
    facet_wrap(~param) + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5)) + 
    theme_bw() + 
    ggtitle("bins") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  
  # Get TSS from bins  ------------------------------------------------------
  
  
  m <- ggplot(params.merge.bins.filt, aes(x = estimate.x, y = estimate.y)) + 
    geom_point(alpha = 0.25) + 
    xlab("logFC norm by total") + 
    ylab("logFC norm by spikeins") + 
    geom_errorbarh(mapping = aes(xmin = estimate.x - se.x, xmax = estimate.x + se.x), alpha = 0.25) + 
    geom_errorbar(mapping = aes(ymin = estimate.y - se.y, ymax = estimate.y + se.y), alpha = 0.25) + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    geom_hline(yintercept = 0, linetype = "dotted") + 
    facet_wrap(~param) + 
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
    theme_bw() + 
    ggtitle("bins") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(m)
  
  # plot with error bars
  for (ctype in ctypes){
    m <- ggplot(params.merge.bins.filt %>% filter(param == ctype & abs(estimate.x) < 5 & abs(estimate.y) < 5), aes(x = estimate.x, y = estimate.y)) + 
      geom_point(alpha = 0.25) + 
      xlab("logFC norm by total") + 
      ylab("logFC norm by spikeins") + 
      geom_errorbarh(mapping = aes(xmin = estimate.x - se.x, xmax = estimate.x + se.x), alpha = 0.25) + 
      geom_errorbar(mapping = aes(ymin = estimate.y - se.y, ymax = estimate.y + se.y), alpha = 0.25) + 
      geom_vline(xintercept = 0, linetype = "dotted") + 
      geom_hline(yintercept = 0, linetype = "dotted") + 
      facet_wrap(~param) + 
      coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5)) + 
      theme_bw() + 
      ggtitle("bins") + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
  }
  
  m <- ggplot(params.merge.bins.filt, aes(x = estimate.x, fill = param)) + 
    geom_density() + 
    facet_wrap(~param) + 
    xlab("logFC norm by spikeins") + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5)) + 
    theme_bw() + 
    ggtitle("bins") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
  m <- ggplot(params.merge.bins.filt, aes(x = estimate.y, fill = param)) + 
    geom_density() + 
    xlab("logFC norm by spikeins") + 
    facet_wrap(~param) + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    coord_cartesian(xlim = c(-5, 5)) + 
    theme_bw() + 
    ggtitle("bins") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom")
  print(m)
  
if (make.plots){
  dev.off()
}

