# Jake Yeung
# Date of Creation: 2020-11-13
# File: ~/projects/scchic/scripts/rstudioserver_analysis/spikeins/BM_merged_with_old/8-differential_expression_analysis_downstream.hiddendomains.R
# description


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
jinf.tss <- "/home/jyeung/hub_oudenaarden/jyeung/data/databases/gene_tss/gene_tss_winsize.50000.bed"
jchromos <- paste("chr", c(seq(19), "X", "Y"), sep = "")


jsettings <- umap.defaults
jsettings$n_neighbors <- 30
jsettings$min_dist <- 0.1
jsettings$random_state <- 123

hubprefix <- "/home/jyeung/hub_oudenaarden"

cbPalette <- c("#696969", "#56B4E9", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#006400",  "#32CD32", "#FFB6C1", "#0b1b7f", "#ff9f7d", "#eb9d01", "#2c2349", "#753187", "#f80597")


# Load fits -------------------------------------------------------------------

# jmark <- "H3K27me3"
# jtype <- "TSS_10000"
jtype <- "hiddendomains"
# jtype <- "TES"
jmark <- "H3K4me3"
jmeth <- "TotalCounts"

jmarks <- c("H3K4me1", "H3K4me3", "H3K27me3"); names(jmarks) <- jmarks


for (jmark in jmarks){
  
  indir <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/pdfs_all/differential_analysis_BM_AllMerged3"
  outf <- file.path(indir, paste0("DE_analysis_", jmark, ".", jtype, ".", jmeth, ".pdf"))
  # if (file.exists(outf)){
  #   print(paste("outf exists, skipping", outf))
  #   next
  # }
  pdf(outf, useDingbats = FALSE)
  
  inf.fits <- paste0("/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/poisson_fits_BM_AllMerged3/poisson_fit_", jtype, ".", jmark, ".2020-11-13.RData")
  assertthat::assert_that(file.exists(inf.fits))
  
  load(inf.fits, v=T)
  
  jrow.names <- names(jfits.lst)
  names(jrow.names) <- jrow.names
  
  # jrow <- grep(jgene, jrow.names, value = TRUE)
  dat.summary <- lapply(jrow.names, function(jrow){
    ctype.effects <- grep("^Cluster", names(jfits.lst[[jrow]]), value = TRUE)
    fits.sub <- unlist(jfits.lst[[jrow]][ctype.effects])
    dat.summary.tmp <- data.frame(param = names(fits.sub), value = fits.sub, rname = jrow, stringsAsFactors = FALSE)
    return(dat.summary.tmp)
  }) %>%
    bind_rows()
  
  
  mall <- ggplot(dat.summary, aes(x = value/ log(2), fill = param)) + 
    geom_density(alpha = 0.25) + 
    scale_fill_manual(values = cbPalette) + 
    theme_bw() + 
    facet_wrap(~param) + 
    geom_vline(xintercept = 0, linetype = "dotted") + 
    ggtitle(jmark, "All Genes") + 
    coord_cartesian(xlim = c(-10, 10)) + 
    xlab("log2FC relative to HSPCs") + 
    theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  print(mall)
  
  # Annotate bins -----------------------------------------------------------
  
  coords.vec <- paste("chr", sapply(unique(dat.summary$rname), function(x) strsplit(x, ";")[[1]][[1]]), sep = "")
  coords.vec.hash <- hash::hash(coords.vec, unique(dat.summary$rname))
  hd.annot <- AnnotateCoordsFromList.GeneWise(coords.vec = coords.vec, inf.tss = jinf.tss, txdb = TxDb.Mmusculus.UCSC.mm10.knownGene, annodb = "org.Mm.eg.db", chromos.keep = jchromos)  
  
  # rename rname 
  hd.annot.summary <- subset(hd.annot$out2.df, select = c(region_coord, tssname, gene)) %>%
    rowwise() %>%
    mutate(rname = AssignHash(x = region_coord, jhash = coords.vec.hash, null.fill = region_coord))
  
  hd.annot.hash <- hash::hash(hd.annot.summary$rname, hd.annot.summary$tssname)
  dat.summary <- dat.summary %>%
    dplyr::rename(rname.orig = rname) %>%
    rowwise() %>%
    mutate(rname = AssignHash(x = rname.orig, jhash = hd.annot.hash, null.fill = rname.orig))
  
  
  
  # Load gene annots --------------------------------------------------------
  
  infrdata <- "/home/jyeung/hub_oudenaarden/jyeung/data/scChiC/from_rstudioserver/rdata_robjs/integrated_analysis_3_marks.setup_for_poisson_regression/integrated_analysis_3_marks.stringentDE.WithHighLowExprs.singlecells.fewerk27me3_TRUE.forPoissonRegression.CountR1only.2020-06-05.smaller.RData"
  assertthat::assert_that(file.exists(infrdata))
  load(infrdata, v=T); rm(tss.mats.filt.fromref.cellfilt)
  
  
  # Add ensembl  ------------------------------------------------------------
  
  dat.summary <- dat.summary %>%
    rowwise() %>%
    mutate(gene = strsplit(rname, ";")[[1]][[2]],
           ens = AssignHash(x = gene, jhash = g2e, null.fill = gene))
  
  jannots <- names(de.ens.sorted.stringent)
  names(jannots) <- jannots
  
  # # jannot <- "Erythroblast"
  # jannot <- "Bcell"
  # jannot <- "HSCs"
  # jannot <- "Neutrophil"
  # jannot <- "HighExprs"
  # jannot <- "LowExprs"
  # 
  
  
  for (jannot in jannots){
    jens <- as.character(de.ens.sorted.stringent[[jannot]])
    m <- ggplot(dat.summary %>% filter(abs(value) < 10) %>% mutate(gene.in.ctype = ens %in% jens), aes(x = value / log(2), fill = gene.in.ctype)) + 
      geom_density(alpha = 0.5) + 
      theme_bw() + 
      facet_wrap(~param) + 
      ggtitle(jmark, paste0("Genes split into two groups by:", jannot, " or not")) + 
      xlab("log2FC relative to HSPCs") + 
      geom_vline(xintercept = 0, linetype = "dotted") + 
      theme(aspect.ratio=1, panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    print(m)
  }
  
  dev.off()
  
  
}





